"""
Маппер идентификаторов метаболитов и реакций на основе MetaNetX.

Использование:
    from metanetx_mapper import MetaNetXMapper

    mapper = MetaNetXMapper("mapping_output/metabolite_mapping.tsv")

    # Получить BiGG id по KEGG
    mapper["C00001"].bigg        # "h2o"
    mapper["C00001"].seed        # "cpd00001"
    mapper["C00001"].mnx         # "MNXM2"
    mapper["C00001"].description # "Water"

    # Все варианты (если ambiguous)
    mapper["C00001"].bigg_all    # ["h2o"]
    mapper["C00341"].bigg_all    # ["grdp", "ipdp"]

    # Проверка наличия
    "C00001" in mapper           # True

    # Итерация
    for kegg_id, entry in mapper.items():
        print(kegg_id, entry.bigg, entry.seed)

    # Стратегии разрешения неоднозначностей
    mapper = MetaNetXMapper("metabolite_mapping.tsv", ambiguity="first")   # первый (дефолт)
    mapper = MetaNetXMapper("metabolite_mapping.tsv", ambiguity="shortest") # самый короткий id
    mapper = MetaNetXMapper("metabolite_mapping.tsv", ambiguity="manual", 
                            overrides={"C00341": {"bigg": "grdp"}})

    # Реакции — точно так же
    rxn_mapper = MetaNetXMapper("mapping_output/reaction_mapping.tsv")
    rxn_mapper["R00200"].bigg    # "PFK"
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class MappingEntry:
    """Одна запись маппинга: KEGG -> MNX -> {bigg, seed, metacyc, ...}"""
    kegg: str
    mnx_id: str
    _bigg_all: list[str] = field(default_factory=list)
    _seed_all: list[str] = field(default_factory=list)
    _metacyc_all: list[str] = field(default_factory=list)
    _rhea_all: list[str] = field(default_factory=list)
    description: str = ""
    ambiguous: str = ""

    # --- Резолверы (устанавливаются маппером) ---
    _resolved_bigg: Optional[str] = field(default=None, repr=False)
    _resolved_seed: Optional[str] = field(default=None, repr=False)
    _resolved_metacyc: Optional[str] = field(default=None, repr=False)
    _resolved_rhea: Optional[str] = field(default=None, repr=False)

    # -- Удобные проперти: один id (resolved) --

    @property
    def bigg(self) -> str:
        return self._resolved_bigg or ""

    @property
    def seed(self) -> str:
        return self._resolved_seed or ""

    @property
    def metacyc(self) -> str:
        return self._resolved_metacyc or ""

    @property
    def rhea(self) -> str:
        return self._resolved_rhea or ""

    # -- Все варианты --

    @property
    def bigg_all(self) -> list[str]:
        return self._bigg_all

    @property
    def seed_all(self) -> list[str]:
        return self._seed_all

    @property
    def metacyc_all(self) -> list[str]:
        return self._metacyc_all

    @property
    def rhea_all(self) -> list[str]:
        return self._rhea_all

    @property
    def is_ambiguous(self) -> bool:
        return bool(self.ambiguous)

    @property
    def mnx(self) -> str:
        return self.mnx_id

    def to_dict(self) -> dict:
        return {
            "kegg": self.kegg,
            "mnx_id": self.mnx_id,
            "bigg": self.bigg,
            "seed": self.seed,
            "metacyc": self.metacyc,
            "rhea": self.rhea,
            "bigg_all": self._bigg_all,
            "seed_all": self._seed_all,
            "description": self.description,
            "ambiguous": self.ambiguous,
        }

    def __str__(self):
        parts = [f"KEGG:{self.kegg}", f"MNX:{self.mnx_id}"]
        if self.bigg:
            parts.append(f"BiGG:{self.bigg}")
        if self.seed:
            parts.append(f"SEED:{self.seed}")
        if self.description:
            parts.append(f'"{self.description}"')
        if self.is_ambiguous:
            parts.append(f"[AMB: {self.ambiguous}]")
        return " | ".join(parts)


# ============================================================
# Стратегии разрешения неоднозначностей
# ============================================================

def _resolve_first(ids: list[str]) -> str:
    """Берём первый."""
    return ids[0] if ids else ""


def _resolve_shortest(ids: list[str]) -> str:
    """Берём самый короткий id (обычно это canonical / human-readable)."""
    if not ids:
        return ""
    return min(ids, key=len)


def _resolve_longest(ids: list[str]) -> str:
    """Берём самый длинный (более специфичный)."""
    if not ids:
        return ""
    return max(ids, key=len)


AMBIGUITY_STRATEGIES = {
    "first": _resolve_first,
    "shortest": _resolve_shortest,
    "longest": _resolve_longest,
}


# ============================================================
# Основной маппер
# ============================================================

class MetaNetXMapper:
    """
    Dict-like маппер KEGG ID -> MappingEntry.

    Параметры:
        filepath:   путь к TSV (metabolite_mapping.tsv или reaction_mapping.tsv)
        ambiguity:  стратегия разрешения неоднозначностей:
                    "first"    — первый id из списка (дефолт)
                    "shortest" — самый короткий id
                    "longest"  — самый длинный id
                    "manual"   — использовать overrides, остальное — first
        overrides:  ручные переопределения для конкретных KEGG ID,
                    формат: {"C00341": {"bigg": "grdp", "seed": "cpd00283"}}
    """

    def __init__(
        self,
        filepath: str,
        ambiguity: str = "first",
        overrides: Optional[dict[str, dict[str, str]]] = None,
    ):
        self._data: dict[str, MappingEntry] = {}
        self._ambiguity = ambiguity
        self._overrides = overrides or {}
        self._filepath = filepath

        if ambiguity not in AMBIGUITY_STRATEGIES and ambiguity != "manual":
            raise ValueError(
                f"Unknown ambiguity strategy: {ambiguity!r}. "
                f"Choose from: {list(AMBIGUITY_STRATEGIES.keys()) + ['manual']}"
            )

        self._resolver = AMBIGUITY_STRATEGIES.get(ambiguity, _resolve_first)
        self._load(filepath)
        self._resolve_all()

    def _parse_pipe_list(self, value: str) -> list[str]:
        """Парсит '|'-separated строку в список, убирает пустые."""
        if not value:
            return []
        return [x.strip() for x in value.split("|") if x.strip()]

    def _load(self, filepath: str):
        """Загружает TSV в _data."""
        with open(filepath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                kegg = row.get("kegg", "").strip()
                if not kegg:
                    continue

                entry = MappingEntry(
                    kegg=kegg,
                    mnx_id=row.get("mnx_id", "").strip(),
                    _bigg_all=self._parse_pipe_list(row.get("bigg", "")),
                    _seed_all=self._parse_pipe_list(row.get("seed", "")),
                    _metacyc_all=self._parse_pipe_list(row.get("metacyc", "")),
                    _rhea_all=self._parse_pipe_list(row.get("rhea", "")),
                    description=row.get("description", "").strip(),
                    ambiguous=row.get("ambiguous", "").strip(),
                )

                # Если KEGG id уже есть — мержим списки (разные MNX могут
                # маппиться на один KEGG, хотя редко)
                if kegg in self._data:
                    existing = self._data[kegg]
                    for attr in ("_bigg_all", "_seed_all", "_metacyc_all", "_rhea_all"):
                        merged = list(dict.fromkeys(
                            getattr(existing, attr) + getattr(entry, attr)
                        ))
                        setattr(existing, attr, merged)
                    if not existing.description and entry.description:
                        existing.description = entry.description
                    if entry.ambiguous:
                        existing.ambiguous = entry.ambiguous
                else:
                    self._data[kegg] = entry

    def _resolve_all(self):
        """Применяет стратегию разрешения ко всем записям."""
        for kegg, entry in self._data.items():
            overrides = self._overrides.get(kegg, {})

            entry._resolved_bigg = (
                overrides.get("bigg") or self._resolver(entry._bigg_all)
            )
            entry._resolved_seed = (
                overrides.get("seed") or self._resolver(entry._seed_all)
            )
            entry._resolved_metacyc = (
                overrides.get("metacyc") or self._resolver(entry._metacyc_all)
            )
            entry._resolved_rhea = (
                overrides.get("rhea") or self._resolver(entry._rhea_all)
            )

    # ---- dict-like interface ----

    def __getitem__(self, kegg_id: str) -> MappingEntry:
        return self._data[kegg_id]

    def __contains__(self, kegg_id: str) -> bool:
        return kegg_id in self._data

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def get(self, kegg_id: str, default=None) -> Optional[MappingEntry]:
        return self._data.get(kegg_id, default)

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    # ---- Утилиты ----

    def ambiguous_entries(self) -> list[MappingEntry]:
        """Возвращает все записи с неоднозначностями."""
        return [e for e in self._data.values() if e.is_ambiguous]

    def missing_bigg(self) -> list[MappingEntry]:
        """KEGG ID без маппинга на BiGG."""
        return [e for e in self._data.values() if not e.bigg]

    def missing_seed(self) -> list[MappingEntry]:
        """KEGG ID без маппинга на SEED."""
        return [e for e in self._data.values() if not e.seed]

    def coverage_stats(self) -> dict:
        """Статистика покрытия."""
        total = len(self._data)
        if total == 0:
            return {"total": 0}
        return {
            "total": total,
            "has_bigg": sum(1 for e in self._data.values() if e.bigg),
            "has_seed": sum(1 for e in self._data.values() if e.seed),
            "has_metacyc": sum(1 for e in self._data.values() if e.metacyc),
            "ambiguous": sum(1 for e in self._data.values() if e.is_ambiguous),
            "bigg_pct": sum(1 for e in self._data.values() if e.bigg) / total * 100,
            "seed_pct": sum(1 for e in self._data.values() if e.seed) / total * 100,
        }

    def reverse_lookup(self, db: str, ext_id: str) -> list[MappingEntry]:
        """
        Обратный поиск: по BiGG/SEED id найти KEGG.
        db: "bigg", "seed", "metacyc", "rhea"
        """
        attr = f"_{db}_all"
        results = []
        for entry in self._data.values():
            if hasattr(entry, attr) and ext_id in getattr(entry, attr):
                results.append(entry)
        return results

    def __repr__(self):
        stats = self.coverage_stats()
        return (
            f"MetaNetXMapper({self._filepath!r}, "
            f"entries={stats.get('total', 0)}, "
            f"bigg={stats.get('bigg_pct', 0):.0f}%, "
            f"seed={stats.get('seed_pct', 0):.0f}%, "
            f"ambiguity={self._ambiguity!r})"
        )


# ============================================================
# Quick demo
# ============================================================

if __name__ == "__main__":
    import sys

    path = sys.argv[1] if len(sys.argv) > 1 else "resources/metabolite_mapping.tsv"
    print(f"Loading {path}...")

    mapper = MetaNetXMapper(path, ambiguity="shortest")
    print(repr(mapper))
    print()

    stats = mapper.coverage_stats()
    for k, v in stats.items():
        print(f"  {k}: {v}")

    print(f"\nАмбигов: {len(mapper.ambiguous_entries())}")
    for e in mapper.ambiguous_entries()[:10]:
        print(f"  {e}")

    print(f"\nПримеры:")
    for kid in list(mapper.keys())[:15]:
        e = mapper[kid]
        print(f"  {e}")