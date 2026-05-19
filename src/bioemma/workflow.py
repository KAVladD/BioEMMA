from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from urllib.request import urlopen

from bioemma.mapper_base import EscherMapper
from bioemma.maps import KeggMap
from bioemma.merger import EscherMerger


@dataclass
class BioEmmaResult:
    """In-memory result plus any files written by the workflow."""

    escher_map: list[dict[str, Any]]
    kegg_reconstruction: dict[str, Any]
    kegg_escher_map: list[dict[str, Any]] | None = None
    fluxes: dict[str, float] | None = None
    paths: dict[str, Path] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)


@dataclass
class BioEmmaBatchResult:
    results: list[BioEmmaResult]
    merged_map: list[dict[str, Any]] | None = None
    paths: dict[str, Path] = field(default_factory=dict)
    summary: dict[str, Any] = field(default_factory=dict)


def normalize_pathway_id(pathway: str | int) -> str:
    value = str(pathway).strip()
    if value.startswith("rn"):
        return value
    return f"rn{value}"


def pathway_url(pathway: str | int) -> str:
    return f"https://rest.kegg.jp/get/{normalize_pathway_id(pathway)}/kgml"


def load_kegg_map(
    *,
    pathway: str | int | None = None,
    kgml: str | Path | None = None,
) -> tuple[KeggMap, dict[str, Any]]:
    if not pathway and not kgml:
        raise ValueError("Either pathway or kgml is required.")

    kegg_map = KeggMap()

    if kgml is not None:
        kgml_path = Path(kgml)
        with kgml_path.open("r", encoding="utf-8") as file:
            kegg_map.read_from_file(file)
        return kegg_map, {"kind": "kgml", "path": str(kgml_path)}

    assert pathway is not None
    pathway_id = normalize_pathway_id(pathway)
    url = pathway_url(pathway_id)
    with urlopen(url) as response:
        kgml_text = response.read().decode("utf-8")
    kegg_map.read_from_file(kgml_text)
    return kegg_map, {"kind": "pathway", "pathway": pathway_id, "url": url}


def load_model(model: Any):
    if hasattr(model, "reactions") and hasattr(model, "metabolites"):
        return model

    _configure_cobra_cache()

    import cobra

    return cobra.io.read_sbml_model(str(model))


def _configure_cobra_cache() -> None:
    cache_dir = os.environ.get("BIOEMMA_COBRA_CACHE_DIR")
    if not cache_dir:
        return

    import appdirs

    original_user_cache_dir = appdirs.user_cache_dir

    def user_cache_dir(appname=None, appauthor=None, *args, **kwargs):
        if appname == "cobrapy" and appauthor == "opencobra":
            return str(Path(cache_dir))
        return original_user_cache_dir(appname, appauthor, *args, **kwargs)

    appdirs.user_cache_dir = user_cache_dir


def summarize_model(model: Any) -> dict[str, Any]:
    """Return stable model metadata for workflow summaries."""

    return {
        "id": getattr(model, "id", None),
        "name": getattr(model, "name", None),
        "reactions": len(getattr(model, "reactions", [])),
        "metabolites": len(getattr(model, "metabolites", [])),
    }


def reconstruct_kegg_map(kegg_map: KeggMap, source: dict[str, Any] | None = None) -> dict[str, Any]:
    metabolites = kegg_map.get_metabolites()
    reactions = kegg_map.get_reactions()
    return {
        "source": source or {},
        "counts": {
            "metabolites": len(metabolites),
            "reactions": len(reactions),
        },
        "metabolites": metabolites,
        "reactions": reactions,
    }


def compute_identifier_coverage(kegg_reconstruction: dict[str, Any]) -> dict[str, Any]:
    namespaces = ("KEGG", "BIGG", "SEED")
    coverage = {}
    for element_type in ("metabolites", "reactions"):
        elements = kegg_reconstruction[element_type]
        total = len(elements)
        coverage[element_type] = {"total": total}
        for namespace in namespaces:
            mapped = sum(
                1
                for element in elements.values()
                if element.get("ids", {}).get(namespace)
            )
            coverage[element_type][namespace] = {
                "mapped": mapped,
                "unmapped": total - mapped,
                "percent": round(mapped / total * 100, 1) if total else 0.0,
            }
    return coverage


def coerce_fluxes(
    model: Any,
    fluxes: Any = None,
    *,
    run_fba: bool = False,
) -> dict[str, float] | None:
    if fluxes is None and run_fba:
        solution = model.optimize()
        fluxes = solution.fluxes

    if fluxes is None:
        return None

    if hasattr(fluxes, "to_dict"):
        fluxes = fluxes.to_dict()

    return {str(key): float(value) for key, value in dict(fluxes).items()}


def build_escher_map(
    model: Any,
    *,
    pathway: str | int | None = None,
    kgml: str | Path | None = None,
    database: str = "BIGG",
    scaling_factor: float = 4,
    axis_epsilon: float = 2,
    remove_orphan_metabolites: bool = False,
    include_kegg_only: bool = False,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    cobra_model = load_model(model)
    kegg_map, source = load_kegg_map(pathway=pathway, kgml=kgml)
    kegg_reconstruction = reconstruct_kegg_map(kegg_map, source)

    mapper = EscherMapper(
        metabolites=kegg_reconstruction["metabolites"],
        reactions=kegg_reconstruction["reactions"],
        database=database,
        scaling_factor=scaling_factor,
        axis_epsilon=axis_epsilon,
        remove_orphan_metabolites=remove_orphan_metabolites,
        include_kegg_only=include_kegg_only,
    )
    escher_map = mapper.build_map(cobra_model)
    kegg_reconstruction["map_stats"] = mapper.map_stats
    return escher_map, kegg_reconstruction


def validate_escher_map(
    escher_map: list[dict[str, Any]],
    *,
    strict_json_keys: bool = False,
) -> dict[str, Any]:
    if not isinstance(escher_map, list) or len(escher_map) != 2:
        raise ValueError("Expected an Escher map shaped as [description, model].")

    description = escher_map[0]
    model = escher_map[1]
    if not isinstance(description, dict) or not isinstance(model, dict):
        raise ValueError("Expected Escher map description and model to be objects.")

    required_description_keys = {
        "map_name",
        "map_id",
        "map_description",
        "homepage",
        "schema",
    }
    required_model_keys = {"nodes", "reactions", "text_labels", "canvas"}
    missing_description_keys = sorted(required_description_keys - set(description.keys()))
    missing_model_keys = sorted(required_model_keys - set(model.keys()))

    nodes = model.get("nodes", {})
    reactions = model.get("reactions", {})
    node_ids = {str(node_id) for node_id in nodes}
    bad_segment_refs = []
    duplicate_segment_ids = []
    seen_segment_ids = set()
    non_string_node_ids = []
    non_string_reaction_ids = []
    non_string_segment_ids = []
    segment_count = 0

    for node_id in nodes:
        if strict_json_keys and not isinstance(node_id, str):
            non_string_node_ids.append(node_id)

    for reaction_id, reaction in reactions.items():
        if strict_json_keys and not isinstance(reaction_id, str):
            non_string_reaction_ids.append(reaction_id)
        for segment_id, segment in reaction.get("segments", {}).items():
            segment_count += 1
            if strict_json_keys and not isinstance(segment_id, str):
                non_string_segment_ids.append(segment_id)
            segment_id_str = str(segment_id)
            if segment_id_str in seen_segment_ids:
                duplicate_segment_ids.append(segment_id_str)
            seen_segment_ids.add(segment_id_str)
            from_node = str(segment.get("from_node_id"))
            to_node = str(segment.get("to_node_id"))
            if from_node not in node_ids or to_node not in node_ids:
                bad_segment_refs.append(
                    {
                        "reaction_id": str(reaction_id),
                        "segment_id": str(segment_id),
                        "from_node_id": from_node,
                        "to_node_id": to_node,
                    }
                )

    return {
        "nodes": len(nodes),
        "reactions": len(reactions),
        "segments": segment_count,
        "missing_description_keys": missing_description_keys,
        "missing_model_keys": missing_model_keys,
        "non_string_node_ids": non_string_node_ids,
        "non_string_reaction_ids": non_string_reaction_ids,
        "non_string_segment_ids": non_string_segment_ids,
        "duplicate_segment_ids": duplicate_segment_ids,
        "bad_segment_refs": bad_segment_refs,
    }


def build_outputs(
    *,
    model: Any,
    pathway: str | int | None = None,
    kgml: str | Path | None = None,
    output_dir: str | Path | None = None,
    map_json_path: str | Path | None = None,
    fluxes: Any = None,
    run_fba: bool = False,
    database: str = "BIGG",
    scaling_factor: float = 4,
    axis_epsilon: float = 2,
    remove_orphan_metabolites: bool = False,
    include_kegg_only: bool = False,
    save_kegg_map: bool = False,
    save_html: bool = False,
    save_png: bool = False,
) -> BioEmmaResult:
    cobra_model = load_model(model)
    kegg_map, source = load_kegg_map(pathway=pathway, kgml=kgml)
    kegg_reconstruction = reconstruct_kegg_map(kegg_map, source)

    mapper = EscherMapper(
        metabolites=kegg_reconstruction["metabolites"],
        reactions=kegg_reconstruction["reactions"],
        database=database,
        scaling_factor=scaling_factor,
        axis_epsilon=axis_epsilon,
        remove_orphan_metabolites=remove_orphan_metabolites,
        include_kegg_only=include_kegg_only,
    )
    escher_map = mapper.build_map(cobra_model)
    kegg_escher_map = None
    if save_kegg_map:
        kegg_mapper = EscherMapper(
            metabolites=kegg_reconstruction["metabolites"],
            reactions=kegg_reconstruction["reactions"],
            database=database,
            scaling_factor=scaling_factor,
            axis_epsilon=axis_epsilon,
        )
        kegg_escher_map = kegg_mapper.build_kegg_map()
    coerced_fluxes = coerce_fluxes(cobra_model, fluxes, run_fba=run_fba)

    paths: dict[str, Path] = {}
    summary = {
        "model": summarize_model(cobra_model),
        "database": database,
        "kegg": kegg_reconstruction["counts"],
        "identifier_coverage": compute_identifier_coverage(kegg_reconstruction),
        "escher": validate_escher_map(escher_map),
        "map_stats": mapper.map_stats,
        "has_fluxes": coerced_fluxes is not None,
    }
    if kegg_escher_map is not None:
        summary["kegg_escher"] = validate_escher_map(kegg_escher_map)

    target_dir: Path | None = None
    if map_json_path is not None:
        paths["escher_map_json"] = Path(map_json_path)
        target_dir = paths["escher_map_json"].parent
    elif output_dir is not None:
        target_dir = Path(output_dir) / _output_slug(pathway=pathway, kgml=kgml)
        paths["escher_map_json"] = target_dir / "escher_map.json"

    if target_dir is not None:
        target_dir.mkdir(parents=True, exist_ok=True)
        paths.setdefault("escher_map_json", target_dir / "escher_map.json")
        paths["kegg_source_reconstruction_json"] = (
            target_dir / "kegg_source_reconstruction.json"
        )
        paths["summary_json"] = target_dir / "summary.json"

        _write_json(paths["escher_map_json"], escher_map)
        _write_json(paths["kegg_source_reconstruction_json"], kegg_reconstruction)
        if kegg_escher_map is not None:
            paths["kegg_escher_map_json"] = target_dir / "kegg_escher_map.json"
            _write_json(paths["kegg_escher_map_json"], kegg_escher_map)

        if coerced_fluxes is not None:
            paths["fluxes_json"] = target_dir / "fluxes.json"
            _write_json(paths["fluxes_json"], coerced_fluxes)

        if save_html:
            paths["escher_map_html"] = target_dir / "escher_map.html"
            _save_html(paths["escher_map_html"], paths["escher_map_json"])
            if kegg_escher_map is not None:
                paths["kegg_escher_map_html"] = target_dir / "kegg_escher_map.html"
                _save_html(paths["kegg_escher_map_html"], paths["kegg_escher_map_json"])
            if coerced_fluxes is not None:
                paths["escher_map_with_fluxes_html"] = (
                    target_dir / "escher_map_with_fluxes.html"
                )
                _save_html(
                    paths["escher_map_with_fluxes_html"],
                    paths["escher_map_json"],
                    model=cobra_model,
                    reaction_data=coerced_fluxes,
                )

        if save_png:
            if "escher_map_html" not in paths:
                paths["escher_map_html"] = target_dir / "escher_map.html"
                _save_html(paths["escher_map_html"], paths["escher_map_json"])
            if kegg_escher_map is not None and "kegg_escher_map_html" not in paths:
                paths["kegg_escher_map_html"] = target_dir / "kegg_escher_map.html"
                _save_html(paths["kegg_escher_map_html"], paths["kegg_escher_map_json"])
            if coerced_fluxes is not None and "escher_map_with_fluxes_html" not in paths:
                paths["escher_map_with_fluxes_html"] = (
                    target_dir / "escher_map_with_fluxes.html"
                )
                _save_html(
                    paths["escher_map_with_fluxes_html"],
                    paths["escher_map_json"],
                    model=cobra_model,
                    reaction_data=coerced_fluxes,
                )
            paths["escher_map_png"] = target_dir / "escher_map.png"
            _save_png(paths["escher_map_html"], paths["escher_map_png"])
            if kegg_escher_map is not None:
                paths["kegg_escher_map_png"] = target_dir / "kegg_escher_map.png"
                _save_png(paths["kegg_escher_map_html"], paths["kegg_escher_map_png"])
            if "escher_map_with_fluxes_html" in paths:
                paths["escher_map_with_fluxes_png"] = (
                    target_dir / "escher_map_with_fluxes.png"
                )
                _save_png(
                    paths["escher_map_with_fluxes_html"],
                    paths["escher_map_with_fluxes_png"],
                )

        summary["paths"] = {key: str(path) for key, path in paths.items()}
        _write_json(paths["summary_json"], summary)

    return BioEmmaResult(
        escher_map=escher_map,
        kegg_reconstruction=kegg_reconstruction,
        kegg_escher_map=kegg_escher_map,
        fluxes=coerced_fluxes,
        paths=paths,
        summary=summary,
    )


def build_many_outputs(
    *,
    model: Any,
    pathways: list[str | int] | None = None,
    kgmls: list[str | Path] | None = None,
    output_dir: str | Path,
    merge: bool = True,
    merged_name: str = "merged_escher_map.json",
    fluxes: Any = None,
    run_fba: bool = False,
    database: str = "BIGG",
    scaling_factor: float = 4,
    axis_epsilon: float = 2,
    remove_orphan_metabolites: bool = False,
    include_kegg_only: bool = False,
    save_kegg_map: bool = False,
    save_html: bool = False,
    save_png: bool = False,
) -> BioEmmaBatchResult:
    if bool(pathways) == bool(kgmls):
        raise ValueError("Provide either pathways or kgmls.")

    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)

    items = list(pathways or kgmls or [])
    results: list[BioEmmaResult] = []
    for item in items:
        kwargs: dict[str, Any]
        if pathways:
            kwargs = {"pathway": item}
        else:
            kwargs = {"kgml": item}

        results.append(
            build_outputs(
                model=model,
                output_dir=output_root,
                fluxes=fluxes,
                run_fba=run_fba,
                database=database,
                scaling_factor=scaling_factor,
                axis_epsilon=axis_epsilon,
                remove_orphan_metabolites=remove_orphan_metabolites,
                include_kegg_only=include_kegg_only,
                save_kegg_map=save_kegg_map,
                save_html=save_html,
                save_png=save_png,
                **kwargs,
            )
        )

    paths: dict[str, Path] = {}
    merged_map = None
    merged_kegg_map = None
    if merge:
        map_paths = [result.paths["escher_map_json"] for result in results]
        merged_map = merge_saved_maps(map_paths)
        paths["merged_escher_map_json"] = output_root / merged_name
        _write_json(paths["merged_escher_map_json"], merged_map)
        if save_html:
            paths["merged_escher_map_html"] = paths[
                "merged_escher_map_json"
            ].with_suffix(".html")
            _save_html(
                paths["merged_escher_map_html"],
                paths["merged_escher_map_json"],
            )
        if save_png:
            if "merged_escher_map_html" not in paths:
                paths["merged_escher_map_html"] = paths[
                    "merged_escher_map_json"
                ].with_suffix(".html")
                _save_html(
                    paths["merged_escher_map_html"],
                    paths["merged_escher_map_json"],
                )
            paths["merged_escher_map_png"] = paths[
                "merged_escher_map_json"
            ].with_suffix(".png")
            _save_png(paths["merged_escher_map_html"], paths["merged_escher_map_png"])
        if save_kegg_map:
            kegg_map_paths = [
                result.paths["kegg_escher_map_json"]
                for result in results
                if "kegg_escher_map_json" in result.paths
            ]
            if kegg_map_paths:
                merged_kegg_map = merge_saved_maps(kegg_map_paths)
                paths["merged_kegg_escher_map_json"] = (
                    output_root / "merged_kegg_escher_map.json"
                )
                _write_json(paths["merged_kegg_escher_map_json"], merged_kegg_map)
                if save_html:
                    paths["merged_kegg_escher_map_html"] = paths[
                        "merged_kegg_escher_map_json"
                    ].with_suffix(".html")
                    _save_html(
                        paths["merged_kegg_escher_map_html"],
                        paths["merged_kegg_escher_map_json"],
                    )
                if save_png:
                    if "merged_kegg_escher_map_html" not in paths:
                        paths["merged_kegg_escher_map_html"] = paths[
                            "merged_kegg_escher_map_json"
                        ].with_suffix(".html")
                        _save_html(
                            paths["merged_kegg_escher_map_html"],
                            paths["merged_kegg_escher_map_json"],
                        )
                    paths["merged_kegg_escher_map_png"] = paths[
                        "merged_kegg_escher_map_json"
                    ].with_suffix(".png")
                    _save_png(
                        paths["merged_kegg_escher_map_html"],
                        paths["merged_kegg_escher_map_png"],
                    )

    summary = {
        "count": len(results),
        "items": [result.summary for result in results],
        "paths": {key: str(path) for key, path in paths.items()},
    }
    if merged_map is not None:
        summary["merged"] = validate_escher_map(merged_map)
    if merged_kegg_map is not None:
        summary["merged_kegg_escher"] = validate_escher_map(merged_kegg_map)

    paths["batch_summary_json"] = output_root / "summary.json"
    _write_json(paths["batch_summary_json"], summary)

    return BioEmmaBatchResult(
        results=results,
        merged_map=merged_map,
        paths=paths,
        summary=summary,
    )


def merge_saved_maps(map_paths: list[str | Path]) -> list[dict[str, Any]]:
    if not map_paths:
        raise ValueError("At least one map path is required.")

    merger = EscherMerger()
    for map_path in map_paths:
        with Path(map_path).open("r", encoding="utf-8") as file:
            merger.append(_compact_escher_map(json.load(file)))
    return merger.merge_maps()


def _compact_escher_map(escher_map: list[dict[str, Any]]) -> list[dict[str, Any]]:
    compacted = json.loads(json.dumps(escher_map))
    model = compacted[1]

    node_id_map = {
        str(old_id): str(new_id)
        for new_id, old_id in enumerate(model.get("nodes", {}).keys())
    }
    model["nodes"] = {
        node_id_map[str(old_id)]: node
        for old_id, node in model.get("nodes", {}).items()
    }

    compacted_reactions = {}
    segment_counter = 0
    for reaction_index, reaction in enumerate(model.get("reactions", {}).values()):
        compacted_segments = {}
        for segment in reaction.get("segments", {}).values():
            segment["from_node_id"] = node_id_map[str(segment["from_node_id"])]
            segment["to_node_id"] = node_id_map[str(segment["to_node_id"])]
            compacted_segments[str(segment_counter)] = segment
            segment_counter += 1
        reaction["segments"] = compacted_segments
        compacted_reactions[str(reaction_index)] = reaction

    model["reactions"] = compacted_reactions
    return compacted


def _output_slug(*, pathway: str | int | None, kgml: str | Path | None) -> str:
    if pathway is not None:
        return normalize_pathway_id(pathway)
    if kgml is not None:
        return Path(kgml).stem
    return "bioemma_map"


def _write_json(path: Path, data: Any) -> None:
    with path.open("w", encoding="utf-8") as file:
        json.dump(data, file, indent=2)


def _save_html(
    path: Path,
    map_json_path: Path,
    *,
    model: Any = None,
    reaction_data: dict[str, float] | None = None,
) -> None:
    try:
        import escher
    except ImportError as exc:
        raise RuntimeError("Saving HTML requires the escher package.") from exc

    kwargs: dict[str, Any] = {"map_json": str(map_json_path)}
    if model is not None:
        kwargs["model"] = model
    if reaction_data is not None:
        kwargs["reaction_data"] = reaction_data

    builder = escher.Builder(**kwargs)
    builder.save_html(str(path))


def _save_png(html_path: Path, png_path: Path) -> None:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as exc:
        raise RuntimeError("Saving PNG requires playwright.") from exc

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1600, "height": 1200})
        page.goto(html_path.resolve().as_uri())
        page.wait_for_load_state("networkidle")
        page.screenshot(path=str(png_path), full_page=True)
        browser.close()
