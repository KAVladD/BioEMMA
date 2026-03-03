#!/usr/bin/env python3
"""
Скрипт для скачивания и обработки кросс-референсных таблиц MetaNetX.
Строит маппинг KEGG <-> BiGG <-> SEED через MNXref namespace
для метаболитов (chem_xref.tsv) и реакций (reac_xref.tsv).

Выход:
  - metabolite_mapping.tsv  (KEGG compound -> MNX -> BiGG -> SEED)
  - reaction_mapping.tsv    (KEGG reaction -> MNX -> BiGG -> SEED)
  - mapping_stats.txt       (статистика покрытия)

Автор: сгенерировано для Escher map pipeline
Дата: 2025
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

# Какие пространства имён нас интересуют
# Формат в chem_xref: "kegg.compound:C00001", "bigg.metabolite:h2o", "seed.compound:cpd00001"
# Формат в reac_xref: "kegg.reaction:R00001", "bigg.reaction:PFK", "seed.reaction:rxn00001"
CHEM_NAMESPACES = {
    "kegg.compound": "kegg",
    "bigg.metabolite": "bigg",
    "seed.compound": "seed",
    "metacyc.compound": "metacyc",  # бонусом, может пригодиться
}

REAC_NAMESPACES = {
    "kegg.reaction": "kegg",
    "bigg.reaction": "bigg",
    "seed.reaction": "seed",
    "metacyc.reaction": "metacyc",
    "rhea": "rhea",
}


# ============================================================
# Скачивание
# ============================================================

def download_file(url: str, dest: str) -> str:
    """Скачивает файл, если ещё не скачан. Возвращает путь."""
    if os.path.exists(dest):
        print(f"  [skip] {dest} уже существует")
        return dest

    print(f"  [download] {url}")
    print(f"             -> {dest}")
    try:
        urllib.request.urlretrieve(url, dest)
        size_mb = os.path.getsize(dest) / (1024 * 1024)
        print(f"  [done] {size_mb:.1f} MB")
    except Exception as e:
        print(f"  [ERROR] Не удалось скачать: {e}")
        sys.exit(1)
    return dest


def download_all():
    """Скачивает все нужные файлы MetaNetX."""
    os.makedirs(DATA_DIR, exist_ok=True)
    paths = {}
    for name, url in FILES.items():
        dest = os.path.join(DATA_DIR, f"{name}.tsv")
        paths[name] = download_file(url, dest)
    return paths


# ============================================================
# Парсинг xref файлов
# ============================================================

def parse_xref(filepath: str, namespaces: dict[str, str]) -> dict[str, dict[str, list[str]]]:
    """
    Парсит chem_xref.tsv или reac_xref.tsv.

    Формат MNXref 4.0:
      chem_xref: XREF \t MNX_ID \t Description
      reac_xref: XREF \t MNX_ID \t Description

    XREF имеет формат "source:id", например "kegg.compound:C00001"

    Возвращает: {mnx_id: {namespace_short: [id1, id2, ...]}}
    Например: {"MNXM1": {"kegg": ["C00080"], "bigg": ["h"], "seed": ["cpd00067"]}}
    """
    mapping = defaultdict(lambda: defaultdict(list))
    skipped = 0
    parsed = 0

    print(f"\n  Парсим {filepath}...")

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            # Пропускаем комменты
            if line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                skipped += 1
                continue

            xref = parts[0]
            mnx_id = parts[1]

            # Парсим "source:id"
            if ":" not in xref:
                skipped += 1
                continue

            # Берём только первое двоеточие (на случай если id содержит ":")
            colon_idx = xref.index(":")
            source = xref[:colon_idx]
            ext_id = xref[colon_idx + 1:]

            if source in namespaces:
                ns_short = namespaces[source]
                mapping[mnx_id][ns_short].append(ext_id)
                parsed += 1
            else:
                skipped += 1

    print(f"  Найдено записей по нужным namespace: {parsed}")
    print(f"  Пропущено (другие namespace/комменты): {skipped}")
    print(f"  Уникальных MNX ID: {len(mapping)}")

    return dict(mapping)


# ============================================================
# Сборка итоговой таблицы
# ============================================================

def build_mapping_table(
    mnx_mapping: dict[str, dict[str, list[str]]],
    primary_ns: str = "kegg",
    target_namespaces: Optional[list[str]] = None,
) -> list[dict]:
    """
    Из маппинга {mnx_id: {ns: [ids]}} строит плоскую таблицу.

    Каждая строка = одна комбинация (primary_id, mnx_id, target_id_1, target_id_2, ...).
    Если для одного MNX есть несколько KEGG / BiGG / SEED — делаем декартово произведение
    (чтобы ничего не потерять), но помечаем неоднозначности.
    """
    if target_namespaces is None:
        target_namespaces = ["bigg", "seed"]

    rows = []

    for mnx_id, ns_ids in mnx_mapping.items():
        # Нужно чтобы был хотя бы primary namespace
        if primary_ns not in ns_ids:
            continue

        primary_ids = ns_ids[primary_ns]

        # Собираем target ids
        target_ids = {}
        for tns in target_namespaces:
            target_ids[tns] = ns_ids.get(tns, [])

        for pid in primary_ids:
            row_base = {
                primary_ns: pid,
                "mnx_id": mnx_id,
            }

            # Если хотя бы один target пустой — всё равно пишем строку
            # Если несколько — раскрываем (но ограничиваем чтобы не взорвалось)
            # Простой подход: join через "|" для множественных значений
            row = dict(row_base)
            for tns in target_namespaces:
                ids = target_ids[tns]
                row[tns] = "|".join(sorted(set(ids))) if ids else ""

            # Пометка неоднозначности
            ambiguous_fields = []
            if len(primary_ids) > 1:
                ambiguous_fields.append(f"{primary_ns}({len(primary_ids)})")
            for tns in target_namespaces:
                if len(target_ids[tns]) > 1:
                    ambiguous_fields.append(f"{tns}({len(target_ids[tns])})")

            row["ambiguous"] = ",".join(ambiguous_fields) if ambiguous_fields else ""

            rows.append(row)

    return rows


def add_descriptions(
    rows: list[dict],
    mnx_mapping: dict[str, dict[str, list[str]]],
    xref_filepath: str,
    primary_ns_prefix: str,
) -> list[dict]:
    """
    Добавляет описания (имена) из xref файла.
    Ищет описание по первому попавшемуся KEGG id для данного MNX.
    """
    # Собираем описания из файла
    descriptions = {}  # xref -> description
    with open(xref_filepath, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                descriptions[parts[0]] = parts[2]

    # Добавляем к строкам
    for row in rows:
        kegg_id = row.get("kegg", "")
        mnx_id = row.get("mnx_id", "")
        desc = ""

        # Пробуем найти описание по kegg id
        kegg_xref = f"{primary_ns_prefix}:{kegg_id}"
        if kegg_xref in descriptions:
            desc = descriptions[kegg_xref]

        # Если не нашли по kegg — попробуем по MNX
        if not desc:
            mnx_xref = f"MNXM:{mnx_id}" if mnx_id.startswith("MNXM") else mnx_id
            # Просто поищем любое описание для этого MNX через bigg
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
# Запись результатов
# ============================================================

def write_tsv(rows: list[dict], filepath: str, columns: list[str]):
    """Записывает список словарей в TSV."""
    os.makedirs(os.path.dirname(filepath) or ".", exist_ok=True)

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        # Сортируем по primary id
        rows_sorted = sorted(rows, key=lambda r: r.get(columns[0], ""))
        for row in rows_sorted:
            writer.writerow(row)

    print(f"  Записано {len(rows)} строк -> {filepath}")


def compute_stats(rows: list[dict], name: str, target_namespaces: list[str]) -> str:
    """Считает статистику покрытия."""
    total = len(rows)
    if total == 0:
        return f"\n=== {name} ===\nНет данных\n"

    lines = [f"\n{'='*60}", f"  {name}", f"{'='*60}"]
    lines.append(f"  Всего строк (KEGG entries): {total}")

    unique_kegg = len(set(r["kegg"] for r in rows))
    lines.append(f"  Уникальных KEGG ID: {unique_kegg}")

    unique_mnx = len(set(r["mnx_id"] for r in rows))
    lines.append(f"  Уникальных MNX ID: {unique_mnx}")

    for tns in target_namespaces:
        has_mapping = sum(1 for r in rows if r.get(tns, ""))
        pct = has_mapping / total * 100
        lines.append(f"  Покрытие {tns}: {has_mapping}/{total} ({pct:.1f}%)")

    ambiguous = sum(1 for r in rows if r.get("ambiguous", ""))
    lines.append(f"  С неоднозначностями: {ambiguous}/{total} ({ambiguous/total*100:.1f}%)")

    # Примеры неоднозначностей
    amb_examples = [r for r in rows if r.get("ambiguous", "")][:5]
    if amb_examples:
        lines.append(f"\n  Примеры неоднозначностей:")
        for r in amb_examples:
            lines.append(f"    KEGG={r['kegg']} MNX={r['mnx_id']} "
                         f"bigg={r.get('bigg','')} seed={r.get('seed','')} "
                         f"[{r['ambiguous']}]")

    # Примеры пробелов (есть KEGG, но нет BiGG или SEED)
    gaps = [r for r in rows if not r.get("bigg", "") or not r.get("seed", "")][:5]
    if gaps:
        lines.append(f"\n  Примеры пробелов:")
        for r in gaps:
            missing = []
            if not r.get("bigg", ""):
                missing.append("bigg")
            if not r.get("seed", ""):
                missing.append("seed")
            lines.append(f"    KEGG={r['kegg']} MNX={r['mnx_id']} "
                         f"missing=[{','.join(missing)}] "
                         f"desc={r.get('description', '')[:50]}")

    return "\n".join(lines)


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 60)
    print("  MetaNetX ID Mapping Builder")
    print("  KEGG <-> MNXref <-> BiGG / SEED")
    print("=" * 60)

    # 1. Скачиваем
    print("\n[1/4] Скачиваем файлы MetaNetX...")
    paths = download_all()

    # 2. Парсим
    print("\n[2/4] Парсим cross-reference файлы...")
    chem_mapping = parse_xref(paths["chem_xref"], CHEM_NAMESPACES)
    reac_mapping = parse_xref(paths["reac_xref"], REAC_NAMESPACES)

    # 3. Собираем таблицы
    print("\n[3/4] Собираем маппинг-таблицы...")

    target_ns = ["bigg", "seed", "metacyc"]

    print("\n  -- Метаболиты --")
    chem_rows = build_mapping_table(chem_mapping, "kegg", target_ns)
    chem_rows = add_descriptions(chem_rows, chem_mapping, paths["chem_xref"],
                                  "kegg.compound")
    print(f"  Итого строк: {len(chem_rows)}")

    print("\n  -- Реакции --")
    reac_target = ["bigg", "seed", "metacyc", "rhea"]
    reac_rows = build_mapping_table(reac_mapping, "kegg", reac_target)
    reac_rows = add_descriptions(reac_rows, reac_mapping, paths["reac_xref"],
                                  "kegg.reaction")
    print(f"  Итого строк: {len(reac_rows)}")

    # 4. Записываем
    print("\n[4/4] Записываем результаты...")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    chem_columns = ["kegg", "mnx_id", "bigg", "seed", "metacyc", "description", "ambiguous"]
    write_tsv(chem_rows,
              os.path.join(OUTPUT_DIR, "metabolite_mapping.tsv"),
              chem_columns)

    reac_columns = ["kegg", "mnx_id", "bigg", "seed", "metacyc", "rhea", "description", "ambiguous"]
    write_tsv(reac_rows,
              os.path.join(OUTPUT_DIR, "reaction_mapping.tsv"),
              reac_columns)

    # Статистика
    stats = []
    stats.append(compute_stats(chem_rows, "METABOLITES (chem_xref)", target_ns))
    stats.append(compute_stats(reac_rows, "REACTIONS (reac_xref)", reac_target))

    stats_text = "\n".join(stats)
    print(stats_text)

    stats_path = os.path.join(OUTPUT_DIR, "mapping_stats.txt")
    with open(stats_path, "w", encoding="utf-8") as f:
        f.write(stats_text)
    print(f"\n  Статистика -> {stats_path}")

    print("\n" + "=" * 60)
    print("  Готово!")
    print(f"  Результаты в папке: {OUTPUT_DIR}/")
    print("=" * 60)


if __name__ == "__main__":
    main()