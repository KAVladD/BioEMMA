#!/usr/bin/env python3
"""
Скрипт для скачивания и обработки кросс-референсных таблиц MetaNetX.
Строит маппинг KEGG <-> BiGG <-> SEED через MNXref namespace
для метаболитов (chem_xref.tsv) и реакций (reac_xref.tsv).

EC-fallback: если для KEGG reaction нет прямого BiGG маппинга через MNX,
ищем BiGG реакцию через общий EC номер (тоже из reac_xref.tsv, без API).

Выход:
  - metabolite_mapping.tsv  (KEGG compound -> MNX -> BiGG -> SEED)
  - reaction_mapping.tsv    (KEGG reaction -> MNX -> EC -> BiGG -> SEED)
  - mapping_stats.txt       (статистика покрытия)
"""

import csv
import os
import sys
import urllib.request
from collections import defaultdict
from typing import Optional

# ============================================================
# Конфиг
# ============================================================

METANETX_BASE = "https://www.metanetx.org/cgi-bin/mnxget/mnxref"
FILES = {
    "chem_xref": f"{METANETX_BASE}/chem_xref.tsv",
    "reac_xref": f"{METANETX_BASE}/reac_xref.tsv",
}

DATA_DIR = "metanetx_data"
OUTPUT_DIR = "resources"

CHEM_NAMESPACES = {
    "kegg.compound": "kegg",
    "bigg.metabolite": "bigg",
    "seed.compound": "seed",
    "metacyc.compound": "metacyc",
}

REAC_NAMESPACES = {
    "kegg.reaction": "kegg",
    "bigg.reaction": "bigg",
    "seed.reaction": "seed",
    "metacyc.reaction": "metacyc",
    "rhea": "rhea",
    "ec-code": "ec",
}


# ============================================================
# Скачивание
# ============================================================

def download_file(url, dest):
    if os.path.exists(dest):
        print(f"  [skip] {dest} already exists")
        return dest
    print(f"  [download] {url}")
    try:
        urllib.request.urlretrieve(url, dest)
        size_mb = os.path.getsize(dest) / (1024 * 1024)
        print(f"  [done] {size_mb:.1f} MB")
    except Exception as e:
        print(f"  [ERROR] Download failed: {e}")
        sys.exit(1)
    return dest


def download_all():
    os.makedirs(DATA_DIR, exist_ok=True)
    paths = {}
    for name, url in FILES.items():
        dest = os.path.join(DATA_DIR, f"{name}.tsv")
        paths[name] = download_file(url, dest)
    return paths


# ============================================================
# Парсинг xref файлов
# ============================================================

def parse_xref(filepath, namespaces):
    """
    Парсит chem_xref.tsv или reac_xref.tsv.
    Возвращает: {mnx_id: {namespace_short: [id1, id2, ...]}}
    """
    mapping = defaultdict(lambda: defaultdict(list))
    skipped = 0
    parsed = 0

    print(f"\n  Parsing {filepath}...")

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                skipped += 1
                continue

            xref = parts[0]
            mnx_id = parts[1]

            if ":" not in xref:
                skipped += 1
                continue

            colon_idx = xref.index(":")
            source = xref[:colon_idx]
            ext_id = xref[colon_idx + 1:]

            if source in namespaces:
                ns_short = namespaces[source]
                mapping[mnx_id][ns_short].append(ext_id)
                parsed += 1
            else:
                skipped += 1

    print(f"  Matched entries: {parsed}")
    print(f"  Skipped: {skipped}")
    print(f"  Unique MNX IDs: {len(mapping)}")

    return dict(mapping)


# ============================================================
# Сборка итоговой таблицы
# ============================================================

def build_mapping_table(mnx_mapping, primary_ns="kegg", target_namespaces=None):
    if target_namespaces is None:
        target_namespaces = ["bigg", "seed"]

    rows = []

    for mnx_id, ns_ids in mnx_mapping.items():
        if primary_ns not in ns_ids:
            continue

        primary_ids = ns_ids[primary_ns]

        target_ids = {}
        for tns in target_namespaces:
            target_ids[tns] = ns_ids.get(tns, [])

        for pid in primary_ids:
            row = {primary_ns: pid, "mnx_id": mnx_id}

            for tns in target_namespaces:
                ids = target_ids[tns]
                row[tns] = "|".join(sorted(set(ids))) if ids else ""

            ambiguous_fields = []
            if len(primary_ids) > 1:
                ambiguous_fields.append(f"{primary_ns}({len(primary_ids)})")
            for tns in target_namespaces:
                if len(target_ids[tns]) > 1:
                    ambiguous_fields.append(f"{tns}({len(target_ids[tns])})")

            row["ambiguous"] = ",".join(ambiguous_fields) if ambiguous_fields else ""
            rows.append(row)

    return rows


def add_descriptions(rows, mnx_mapping, xref_filepath, primary_ns_prefix):
    descriptions = {}
    with open(xref_filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                descriptions[parts[0]] = parts[2]

    for row in rows:
        kegg_id = row.get("kegg", "")
        mnx_id = row.get("mnx_id", "")
        desc = ""

        kegg_xref = f"{primary_ns_prefix}:{kegg_id}"
        if kegg_xref in descriptions:
            desc = descriptions[kegg_xref]

        if not desc:
            bigg_ids = row.get("bigg", "").split("|")
            for bid in bigg_ids:
                if bid:
                    bxref = f"bigg.metabolite:{bid}" if "MNXM" in mnx_id else f"bigg.reaction:{bid}"
                    if bxref in descriptions:
                        desc = descriptions[bxref]
                        break

        row["description"] = desc

    return rows


# ============================================================
# EC-based fallback (offline, из reac_xref.tsv)
# ============================================================

def build_ec_fallback_offline(mnx_mapping, reac_rows):
    """
    Для реакций без BiGG: ищет BiGG через общий EC номер.
    
    Проблема: EC и KEGG reaction часто лежат на РАЗНЫХ MNX ID.
    Поэтому строим два индекса:
      1. kegg_rxn -> set(EC) через ВСЕ MNX (не только тот что в строке)
      2. EC -> set(BiGG) через ВСЕ MNX
    Потом: строка без BiGG → её KEGG → EC из индекса 1 → BiGG из индекса 2.
    
    Также обогащает колонку 'ec' у всех строк.
    """
    # Индекс 1: kegg_reaction -> set(EC)
    # Путь: kegg.reaction -> MNXR_A, ec-code -> MNXR_A (если повезло),
    # но чаще: kegg -> MNXR_A, ec -> MNXR_B. Связь через то что оба
    # указывают на "ту же реакцию". В MetaNetX это решается через
    # deprecated/merged MNX IDs. Но в xref они могут быть на разных MNX.
    #
    # Простой путь: собираем kegg->ec НЕ через MNX напрямую,
    # а через любой общий путь. Для этого:
    # - Для каждого MNX собираем все kegg и все ec
    # - Если MNX имеет и kegg и ec — связываем
    kegg_to_ec = defaultdict(set)
    ec_to_bigg = defaultdict(set)

    for mnx_id, ns_ids in mnx_mapping.items():
        keggs = ns_ids.get("kegg", [])
        ecs = ns_ids.get("ec", [])
        biggs = ns_ids.get("bigg", [])

        # kegg -> ec (через общий MNX)
        for k in keggs:
            for ec in ecs:
                kegg_to_ec[k].add(ec)

        # ec -> bigg (через общий MNX)
        for ec in ecs:
            for b in biggs:
                ec_to_bigg[ec].add(b)

    print(f"  KEGG->EC index: {len(kegg_to_ec)} KEGG reactions with EC")
    print(f"  EC->BiGG index: {len(ec_to_bigg)} EC numbers with BiGG")

    # Обогащаем ВСЕ строки EC (даже если уже есть bigg)
    ec_enriched = 0
    for row in reac_rows:
        kegg_id = row.get("kegg", "")
        # Мержим EC из строки (от build_mapping_table) и из индекса
        existing_ec = set(e.strip() for e in row.get("ec", "").split("|") if e.strip())
        from_index = kegg_to_ec.get(kegg_id, set())
        all_ec = existing_ec | from_index
        if all_ec:
            row["ec"] = "|".join(sorted(all_ec))
            if from_index - existing_ec:
                ec_enriched += 1

    if ec_enriched:
        print(f"  EC enriched from index: {ec_enriched} rows got extra EC numbers")

    # Фоллбэк BiGG
    fallback_count = 0
    for row in reac_rows:
        if row.get("bigg", ""):
            continue

        ec_str = row.get("ec", "")
        if not ec_str:
            continue

        ecs = [e.strip() for e in ec_str.split("|") if e.strip()]
        fallback_bigg = set()
        for ec in ecs:
            fallback_bigg.update(ec_to_bigg.get(ec, set()))

        if fallback_bigg:
            row["bigg"] = "|".join(sorted(fallback_bigg))
            amb = row.get("ambiguous", "")
            ec_note = "ec_fallback"
            if len(fallback_bigg) > 1:
                ec_note += f"(bigg:{len(fallback_bigg)})"
            row["ambiguous"] = f"{amb},{ec_note}" if amb else ec_note
            fallback_count += 1

    print(f"  EC fallback: {fallback_count} reactions got BiGG via EC")
    return reac_rows


# ============================================================
# Запись результатов
# ============================================================

def write_tsv(rows, filepath, columns):
    os.makedirs(os.path.dirname(filepath) or ".", exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        rows_sorted = sorted(rows, key=lambda r: r.get(columns[0], ""))
        for row in rows_sorted:
            writer.writerow(row)
    print(f"  Written {len(rows)} rows -> {filepath}")


def compute_stats(rows, name, target_namespaces):
    total = len(rows)
    if total == 0:
        return f"\n=== {name} ===\nNo data\n"

    lines = [f"\n{'='*60}", f"  {name}", f"{'='*60}"]
    lines.append(f"  Total rows (KEGG entries): {total}")
    lines.append(f"  Unique KEGG IDs: {len(set(r['kegg'] for r in rows))}")
    lines.append(f"  Unique MNX IDs: {len(set(r['mnx_id'] for r in rows))}")

    for tns in target_namespaces:
        has_mapping = sum(1 for r in rows if r.get(tns, ""))
        pct = has_mapping / total * 100
        lines.append(f"  Coverage {tns}: {has_mapping}/{total} ({pct:.1f}%)")

    ambiguous = sum(1 for r in rows if r.get("ambiguous", ""))
    lines.append(f"  Ambiguous: {ambiguous}/{total} ({ambiguous/total*100:.1f}%)")

    ec_fb = sum(1 for r in rows if "ec_fallback" in r.get("ambiguous", ""))
    if ec_fb:
        lines.append(f"  EC fallback: {ec_fb}")

    amb_examples = [r for r in rows if r.get("ambiguous", "")][:5]
    if amb_examples:
        lines.append(f"\n  Ambiguity examples:")
        for r in amb_examples:
            lines.append(f"    KEGG={r['kegg']} MNX={r['mnx_id']} "
                         f"bigg={r.get('bigg','')} [{r['ambiguous']}]")

    gaps = [r for r in rows if not r.get("bigg", "") or not r.get("seed", "")][:5]
    if gaps:
        lines.append(f"\n  Gap examples:")
        for r in gaps:
            missing = []
            if not r.get("bigg", ""):
                missing.append("bigg")
            if not r.get("seed", ""):
                missing.append("seed")
            lines.append(f"    KEGG={r['kegg']} MNX={r['mnx_id']} "
                         f"missing=[{','.join(missing)}]")

    return "\n".join(lines)


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("  MetaNetX ID Mapping Builder")
    print("  KEGG <-> MNXref <-> BiGG / SEED + EC fallback")
    print("=" * 60)

    # 1. Download
    print("\n[1/4] Downloading MetaNetX files...")
    paths = download_all()

    # 2. Parse
    print("\n[2/4] Parsing cross-reference files...")
    chem_mapping = parse_xref(paths["chem_xref"], CHEM_NAMESPACES)
    reac_mapping = parse_xref(paths["reac_xref"], REAC_NAMESPACES)

    # 3. Build tables
    print("\n[3/4] Building mapping tables...")

    target_ns = ["bigg", "seed", "metacyc"]

    print("\n  -- Metabolites --")
    chem_rows = build_mapping_table(chem_mapping, "kegg", target_ns)
    chem_rows = add_descriptions(chem_rows, chem_mapping, paths["chem_xref"],
                                  "kegg.compound")
    print(f"  Total rows: {len(chem_rows)}")

    print("\n  -- Reactions --")
    reac_target = ["bigg", "seed", "metacyc", "rhea", "ec"]
    reac_rows = build_mapping_table(reac_mapping, "kegg", reac_target)
    reac_rows = add_descriptions(reac_rows, reac_mapping, paths["reac_xref"],
                                  "kegg.reaction")
    print(f"  Total rows: {len(reac_rows)}")

    # 3.5. EC fallback (offline)
    print("\n  -- EC fallback (offline) --")
    reac_rows = build_ec_fallback_offline(reac_mapping, reac_rows)

    # 4. Write
    print("\n[4/4] Writing results...")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    chem_columns = ["kegg", "mnx_id", "bigg", "seed", "metacyc",
                    "description", "ambiguous"]
    write_tsv(chem_rows,
              os.path.join(OUTPUT_DIR, "metabolite_mapping.tsv"),
              chem_columns)

    reac_columns = ["kegg", "mnx_id", "ec", "bigg", "seed", "metacyc",
                    "rhea", "description", "ambiguous"]
    write_tsv(reac_rows,
              os.path.join(OUTPUT_DIR, "reaction_mapping.tsv"),
              reac_columns)

    # Stats
    stats = []
    stats.append(compute_stats(chem_rows, "METABOLITES (chem_xref)", target_ns))
    stats.append(compute_stats(reac_rows, "REACTIONS (reac_xref)", reac_target))

    stats_text = "\n".join(stats)
    print(stats_text)

    stats_path = os.path.join(OUTPUT_DIR, "mapping_stats.txt")
    with open(stats_path, "w", encoding="utf-8") as f:
        f.write(stats_text)
    print(f"\n  Stats -> {stats_path}")

    print("\n" + "=" * 60)
    print("  Done!")
    print(f"  Results in: {OUTPUT_DIR}/")
    print("=" * 60)


if __name__ == "__main__":
    main()