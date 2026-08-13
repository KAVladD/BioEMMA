#!/usr/bin/env python3
"""
Скрипт для скачивания и обработки кросс-референсных таблиц MetaNetX.
Строит маппинг KEGG <-> BiGG <-> SEED через MNXref namespace
для метаболитов (chem_xref.tsv) и реакций (reac_xref.tsv).

EC-fallback: если для KEGG reaction нет прямого BiGG/SEED маппинга через MNX,
ищем BiGG или SEED реакцию через общий EC номер и проверяем пересечение
метаболитов-участников реакции (тоже из reac_xref.tsv/reac_prop.tsv, без API).

Выход:
  - metabolite_mapping.tsv  (KEGG compound -> MNX -> BiGG -> SEED)
  - reaction_mapping.tsv    (KEGG reaction -> MNX -> EC -> BiGG -> SEED)
  - mapping_stats.txt       (статистика покрытия)
"""

import csv
import os
import re
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
    "reac_prop": f"{METANETX_BASE}/reac_prop.tsv",
}

DATA_DIR = os.environ.get("BIOEMMA_METANETX_DATA_DIR", "metanetx_data")
OUTPUT_DIR = os.environ.get("BIOEMMA_MAPPING_OUTPUT_DIR", "resources")
EC_FALLBACK_MIN_PARTICIPANT_OVERLAP = 0.5

REACTION_PARTICIPANT_RE = re.compile(r"([A-Za-z0-9_.-]+)@MNXD[0-9A-Za-z_]+")

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


def parse_reaction_participants(equation):
    return sorted(set(REACTION_PARTICIPANT_RE.findall(equation)))


def split_pipe(value):
    return [part.strip() for part in value.split("|") if part.strip()]


def has_ambiguous_marker(ambiguous, marker):
    return re.search(rf"(?:^|,)\s*{re.escape(marker)}\(", ambiguous) is not None


def parse_reaction_props(filepath):
    """Return reaction EC numbers and participant sets from MetaNetX reac_prop.tsv."""
    mnx_to_props = defaultdict(lambda: {"ec": [], "participants": []})
    parsed = 0

    print(f"\n  Parsing reaction properties from {filepath}...")

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) < 4:
                continue

            mnx_id = parts[0]
            participants = parse_reaction_participants(parts[1])
            if participants:
                mnx_to_props[mnx_id]["participants"] = participants

            ec_field = parts[3].strip()
            if not ec_field:
                continue

            for ec in ec_field.replace("|", ";").split(";"):
                ec = ec.strip()
                if ec and ec not in mnx_to_props[mnx_id]["ec"]:
                    mnx_to_props[mnx_id]["ec"].append(ec)
                    parsed += 1

    ec_reactions = sum(1 for props in mnx_to_props.values() if props["ec"])
    participant_reactions = sum(
        1 for props in mnx_to_props.values() if props["participants"]
    )
    print(f"  MNX reactions with participants: {participant_reactions}")
    print(f"  MNX reactions with EC: {ec_reactions}")
    print(f"  EC assignments: {parsed}")
    return dict(mnx_to_props)


def add_reaction_props_from_props(reac_mapping, reac_prop_filepath):
    """Merge EC numbers and participants from reac_prop.tsv into reaction xrefs."""
    mnx_to_props = parse_reaction_props(reac_prop_filepath)
    enriched_mnx = 0
    enriched_ec = 0
    enriched_participants = 0

    for mnx_id, props in mnx_to_props.items():
        if mnx_id not in reac_mapping:
            continue

        participants = props.get("participants", [])
        if participants:
            reac_mapping[mnx_id]["participants"] = participants
            enriched_participants += 1

        ecs = props.get("ec", [])
        existing = set(reac_mapping[mnx_id].get("ec", []))
        new_ecs = [ec for ec in ecs if ec not in existing]
        if not new_ecs:
            continue

        reac_mapping[mnx_id]["ec"].extend(new_ecs)
        enriched_mnx += 1
        enriched_ec += len(new_ecs)

    print(f"  Reaction MNX entries enriched with EC: {enriched_mnx}")
    print(f"  EC values added to xref mapping: {enriched_ec}")
    print(f"  Reaction MNX entries enriched with participants: {enriched_participants}")
    return reac_mapping


# ============================================================
# EC-based fallback (offline, из reac_xref.tsv)
# ============================================================

def participant_overlap(source_participants, candidate_participants):
    if not source_participants or not candidate_participants:
        return 0.0
    return len(source_participants & candidate_participants) / len(source_participants)


def build_ec_fallback_offline(
    mnx_mapping,
    reac_rows,
    min_participant_overlap=EC_FALLBACK_MIN_PARTICIPANT_OVERLAP,
):
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
    ec_to_bigg = defaultdict(list)

    for mnx_id, ns_ids in mnx_mapping.items():
        keggs = ns_ids.get("kegg", [])
        ecs = ns_ids.get("ec", [])
        biggs = ns_ids.get("bigg", [])
        participants = set(ns_ids.get("participants", []))

        # kegg -> ec (через общий MNX)
        for k in keggs:
            for ec in ecs:
                kegg_to_ec[k].add(ec)

        # ec -> bigg (через общий MNX)
        for ec in ecs:
            for b in biggs:
                ec_to_bigg[ec].append(
                    {
                        "bigg": b,
                        "mnx_id": mnx_id,
                        "participants": participants,
                    }
                )

    print(f"  KEGG->EC index: {len(kegg_to_ec)} KEGG reactions with EC")
    print(f"  EC->BiGG index: {len(ec_to_bigg)} EC numbers with BiGG")
    print(
        "  EC fallback participant overlap threshold: "
        f"{min_participant_overlap:.0%}"
    )

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
    rejected_by_participants = 0
    discarded_by_best_overlap = 0
    missing_participants = 0
    for row in reac_rows:
        if row.get("bigg", ""):
            continue

        ec_str = row.get("ec", "")
        if not ec_str:
            continue

        source_participants = set(
            mnx_mapping.get(row.get("mnx_id", ""), {}).get("participants", [])
        )
        if not source_participants:
            missing_participants += 1
            continue

        ecs = [e.strip() for e in ec_str.split("|") if e.strip()]
        fallback_scores = {}
        for ec in ecs:
            for candidate in ec_to_bigg.get(ec, []):
                overlap = participant_overlap(
                    source_participants,
                    candidate["participants"],
                )
                if overlap < min_participant_overlap:
                    rejected_by_participants += 1
                    continue

                bigg_id = candidate["bigg"]
                fallback_scores[bigg_id] = max(
                    fallback_scores.get(bigg_id, 0.0),
                    overlap,
                )

        if fallback_scores:
            best_overlap = max(fallback_scores.values())
            fallback_bigg = sorted(
                bigg_id
                for bigg_id, score in fallback_scores.items()
                if score == best_overlap
            )
            discarded_by_best_overlap += len(fallback_scores) - len(fallback_bigg)
            row["bigg"] = "|".join(fallback_bigg)
            amb = row.get("ambiguous", "")
            ec_note = f"ec_fallback(participants>={min_participant_overlap:.2f}"
            if len(fallback_bigg) > 1:
                ec_note += f",bigg:{len(fallback_bigg)}"
            ec_note += ")"
            row["ambiguous"] = f"{amb},{ec_note}" if amb else ec_note
            fallback_count += 1

    print(f"  BiGG EC fallback: {fallback_count} reactions got BiGG via EC")
    print(
        "  BiGG EC fallback candidates rejected by participant overlap: "
        f"{rejected_by_participants}"
    )
    if discarded_by_best_overlap:
        print(
            "  BiGG EC fallback lower-overlap BiGG IDs discarded: "
            f"{discarded_by_best_overlap}"
        )
    if missing_participants:
        print(
            "  BiGG EC fallback rows skipped without participant data: "
            f"{missing_participants}"
        )
    return reac_rows


def build_seed_ec_fallback_offline(
    mnx_mapping,
    reac_rows,
    min_participant_overlap=EC_FALLBACK_MIN_PARTICIPANT_OVERLAP,
):
    """
    Для реакций без SEED: ищет SEED через общий EC номер.

    Кандидаты берутся по всем MNX-реакциям с тем же EC. Кандидат проходит
    только если его набор участников покрывает исходную KEGG-реакцию минимум
    на min_participant_overlap. Если кандидатов несколько, сохраняются только
    кандидаты с максимальным overlap.
    """
    ec_to_seed = defaultdict(list)
    seen = set()

    for mnx_id, ns_ids in mnx_mapping.items():
        seeds = ns_ids.get("seed", [])
        ecs = ns_ids.get("ec", [])
        if not seeds or not ecs:
            continue

        participants = set(ns_ids.get("participants", []))
        for ec in ecs:
            for seed_id in seeds:
                key = (ec, seed_id, mnx_id)
                if key in seen:
                    continue
                seen.add(key)
                ec_to_seed[ec].append(
                    {
                        "seed": seed_id,
                        "mnx_id": mnx_id,
                        "participants": participants,
                    }
                )

    print(f"  EC->SEED index: {len(ec_to_seed)} EC numbers with SEED")

    fallback_count = 0
    rejected_by_participants = 0
    discarded_by_best_overlap = 0
    missing_participants = 0
    no_ec = 0

    for row in reac_rows:
        if row.get("seed", ""):
            continue

        ecs = split_pipe(row.get("ec", ""))
        if not ecs:
            no_ec += 1
            continue

        source_participants = set(
            mnx_mapping.get(row.get("mnx_id", ""), {}).get("participants", [])
        )
        if not source_participants:
            missing_participants += 1
            continue

        fallback_scores = {}
        for ec in ecs:
            for candidate in ec_to_seed.get(ec, []):
                overlap = participant_overlap(
                    source_participants,
                    candidate["participants"],
                )
                if overlap < min_participant_overlap:
                    rejected_by_participants += 1
                    continue

                seed_id = candidate["seed"]
                fallback_scores[seed_id] = max(
                    fallback_scores.get(seed_id, 0.0),
                    overlap,
                )

        if not fallback_scores:
            continue

        best_overlap = max(fallback_scores.values())
        fallback_seed = sorted(
            seed_id
            for seed_id, score in fallback_scores.items()
            if score == best_overlap
        )
        discarded_by_best_overlap += len(fallback_scores) - len(fallback_seed)
        row["seed"] = "|".join(fallback_seed)

        amb = row.get("ambiguous", "")
        seed_note = f"seed_ec_fallback(participants>={min_participant_overlap:.2f}"
        if len(fallback_seed) > 1:
            seed_note += f",seed:{len(fallback_seed)}"
        seed_note += ")"
        row["ambiguous"] = f"{amb},{seed_note}" if amb else seed_note
        fallback_count += 1

    print(f"  SEED EC fallback: {fallback_count} reactions got SEED via EC")
    print(
        "  SEED EC fallback candidates rejected by participant overlap: "
        f"{rejected_by_participants}"
    )
    if discarded_by_best_overlap:
        print(
            "  SEED EC fallback lower-overlap SEED IDs discarded: "
            f"{discarded_by_best_overlap}"
        )
    if missing_participants:
        print(
            "  SEED EC fallback rows skipped without participant data: "
            f"{missing_participants}"
        )
    if no_ec:
        print(f"  SEED EC fallback rows skipped without EC: {no_ec}")
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

    bigg_ec_fb = sum(
        1 for r in rows if has_ambiguous_marker(r.get("ambiguous", ""), "ec_fallback")
    )
    seed_ec_fb = sum(
        1
        for r in rows
        if has_ambiguous_marker(r.get("ambiguous", ""), "seed_ec_fallback")
    )
    if bigg_ec_fb:
        lines.append(f"  BiGG EC fallback: {bigg_ec_fb}")
    if seed_ec_fb:
        lines.append(f"  SEED EC fallback: {seed_ec_fb}")

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
    print("  KEGG <-> MNXref <-> BiGG / SEED + EC fallbacks")
    print("=" * 60)

    # 1. Download
    print("\n[1/4] Downloading MetaNetX files...")
    paths = download_all()

    # 2. Parse
    print("\n[2/4] Parsing cross-reference files...")
    chem_mapping = parse_xref(paths["chem_xref"], CHEM_NAMESPACES)
    reac_mapping = parse_xref(paths["reac_xref"], REAC_NAMESPACES)
    reac_mapping = add_reaction_props_from_props(reac_mapping, paths["reac_prop"])

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

    # 3.5. EC fallbacks (offline)
    print("\n  -- BiGG EC fallback (offline) --")
    reac_rows = build_ec_fallback_offline(reac_mapping, reac_rows)
    print("\n  -- SEED EC fallback (offline) --")
    reac_rows = build_seed_ec_fallback_offline(reac_mapping, reac_rows)

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
    stats.append(compute_stats(reac_rows, "REACTIONS (reac_xref + reac_prop)", reac_target))

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
