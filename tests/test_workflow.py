import json
import os
import subprocess
import sys
from pathlib import Path

from bioemma.workflow import build_many_outputs, build_outputs, validate_escher_map


ROOT = Path(__file__).resolve().parents[1]
KGML = ROOT / "tests" / "data" / "kgml" / "rn00010.xml"
KGML_TCA = ROOT / "tests" / "data" / "kgml" / "rn00020.xml"
MODEL = ROOT / "tests" / "data" / "models" / "e_coli_core.xml"


def _set_cobra_cache(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("LOCALAPPDATA", str(cache_root))
    monkeypatch.setenv("APPDATA", str(cache_root))
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache_root))
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))


def _assert_valid_saved_escher_map(path):
    escher_map = json.loads(path.read_text(encoding="utf-8"))
    validation = validate_escher_map(escher_map, strict_json_keys=True)
    assert validation["missing_description_keys"] == []
    assert validation["missing_model_keys"] == []
    assert validation["non_string_node_ids"] == []
    assert validation["non_string_reaction_ids"] == []
    assert validation["non_string_segment_ids"] == []
    assert validation["duplicate_segment_ids"] == []
    assert validation["bad_segment_refs"] == []
    return escher_map


def test_build_outputs_returns_and_saves_core_artifacts(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        output_dir=tmp_path / "out",
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
        save_kegg_map=True,
    )

    assert result.escher_map[0]["schema"] == "https://escher.github.io/escher/jsonschema/1-0-0#"
    assert result.kegg_reconstruction["counts"]["metabolites"] > 0
    assert result.kegg_reconstruction["counts"]["reactions"] > 0

    validation = validate_escher_map(result.escher_map)
    assert validation["nodes"] > 0
    assert validation["reactions"] < result.kegg_reconstruction["counts"]["reactions"]
    assert validation["segments"] > 0
    assert validation["bad_segment_refs"] == []

    assert result.paths["escher_map_json"].is_file()
    assert result.paths["kegg_escher_map_json"].is_file()
    assert result.paths["kegg_source_reconstruction_json"].is_file()
    assert result.paths["summary_json"].is_file()

    saved_map = _assert_valid_saved_escher_map(result.paths["escher_map_json"])
    saved_kegg_map = _assert_valid_saved_escher_map(result.paths["kegg_escher_map_json"])
    saved_reconstruction = json.loads(
        result.paths["kegg_source_reconstruction_json"].read_text(encoding="utf-8")
    )
    saved_summary = json.loads(result.paths["summary_json"].read_text(encoding="utf-8"))
    assert validate_escher_map(saved_map)["nodes"] == validation["nodes"]
    assert validate_escher_map(saved_kegg_map)["reactions"] == result.kegg_reconstruction[
        "counts"
    ]["reactions"]
    assert saved_reconstruction["counts"] == result.kegg_reconstruction["counts"]
    assert saved_summary["map_stats"] == result.summary["map_stats"]
    assert saved_summary["kegg_escher"]["bad_segment_refs"] == []
    assert saved_summary["model"]["id"] == "e_coli_core"
    assert saved_summary["model"]["reactions"] == 95
    assert saved_summary["model"]["metabolites"] == 72
    assert saved_summary["database"] == "BIGG"
    assert saved_summary["identifier_coverage"]["reactions"]["BIGG"]["mapped"] > 0
    assert saved_summary["identifier_coverage"]["metabolites"]["SEED"]["mapped"] > 0

    stages = {stage["name"]: stage for stage in result.summary["map_stats"]["stages"]}
    assert "kegg_layout" in stages
    assert "model_reaction_filter" in stages
    assert "model_metabolite_filter" in stages
    assert "secondary_metabolite_addition" in stages
    assert "final_layout" in stages
    assert stages["model_reaction_filter"]["change"]["reactions"]["removed"] > 0
    assert stages["secondary_metabolite_addition"]["change"]["nodes"]["added"] > 0


def test_build_outputs_can_include_kegg_only_elements(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
        include_kegg_only=True,
    )

    validation = validate_escher_map(result.escher_map)
    assert validation["reactions"] == result.kegg_reconstruction["counts"]["reactions"]
    assert validation["bad_segment_refs"] == []


def test_cli_build_writes_workflow_outputs(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "bioemma.cli",
            "build",
            "--model",
            str(MODEL),
            "--kgml",
            str(KGML),
            "--output-dir",
            str(tmp_path / "cli-out"),
            "--database",
            "BIGG",
            "--scaling-factor",
            "5",
            "--axis-epsilon",
            "10",
            "--save-kegg-map",
            "--map-stats",
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )

    output_dir = tmp_path / "cli-out" / "rn00010"
    assert (output_dir / "escher_map.json").is_file()
    assert (output_dir / "kegg_escher_map.json").is_file()
    assert (output_dir / "kegg_source_reconstruction.json").is_file()
    assert (output_dir / "summary.json").is_file()
    assert "escher_map_json:" in completed.stdout
    assert "kegg_escher_map_json:" in completed.stdout
    assert "map_stats:" in completed.stdout
    assert "model_reaction_filter:" in completed.stdout
    assert "+0/-" in completed.stdout


def test_build_many_outputs_merges_saved_maps(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    result = build_many_outputs(
        model=MODEL,
        kgmls=[KGML, KGML_TCA],
        output_dir=tmp_path / "batch-out",
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
    )

    assert len(result.results) == 2
    assert (tmp_path / "batch-out" / "rn00010" / "escher_map.json").is_file()
    assert (tmp_path / "batch-out" / "rn00020" / "escher_map.json").is_file()
    assert result.paths["merged_escher_map_json"].is_file()

    validation = validate_escher_map(result.merged_map)
    assert validation["nodes"] > 0
    assert validation["reactions"] > 0
    assert validation["bad_segment_refs"] == []
    assert validation["duplicate_segment_ids"] == []
    _assert_valid_saved_escher_map(result.paths["merged_escher_map_json"])


def test_build_many_outputs_merges_model_and_kegg_maps_separately(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    result = build_many_outputs(
        model=MODEL,
        kgmls=[KGML, KGML_TCA],
        output_dir=tmp_path / "batch-kegg-out",
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
        save_kegg_map=True,
    )

    model_merged = _assert_valid_saved_escher_map(result.paths["merged_escher_map_json"])
    kegg_merged = _assert_valid_saved_escher_map(
        result.paths["merged_kegg_escher_map_json"]
    )
    assert len(model_merged[1]["reactions"]) < len(kegg_merged[1]["reactions"])
    assert result.summary["merged"]["reactions"] < result.summary["merged_kegg_escher"][
        "reactions"
    ]


def test_build_many_outputs_can_render_merged_maps(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    def fake_save_html(path, _map_json_path, **_kwargs):
        path.write_text("<html></html>", encoding="utf-8")

    def fake_save_png(_html_path, png_path):
        png_path.write_bytes(b"png")

    monkeypatch.setattr("bioemma.workflow._save_html", fake_save_html)
    monkeypatch.setattr("bioemma.workflow._save_png", fake_save_png)

    result = build_many_outputs(
        model=MODEL,
        kgmls=[KGML, KGML_TCA],
        output_dir=tmp_path / "batch-render-out",
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
        save_kegg_map=True,
        save_html=True,
        save_png=True,
    )

    for key in [
        "merged_escher_map_html",
        "merged_escher_map_png",
        "merged_kegg_escher_map_html",
        "merged_kegg_escher_map_png",
    ]:
        assert result.paths[key].is_file()


def test_cli_build_accepts_multiple_kgmls_and_merges(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "bioemma.cli",
            "build",
            "--model",
            str(MODEL),
            "--kgml",
            str(KGML),
            str(KGML_TCA),
            "--output-dir",
            str(tmp_path / "cli-batch-out"),
            "--database",
            "BIGG",
            "--scaling-factor",
            "5",
            "--axis-epsilon",
            "10",
            "--save-kegg-map",
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )

    output_dir = tmp_path / "cli-batch-out"
    assert (output_dir / "rn00010" / "escher_map.json").is_file()
    assert (output_dir / "rn00010" / "kegg_escher_map.json").is_file()
    assert (output_dir / "rn00020" / "escher_map.json").is_file()
    assert (output_dir / "rn00020" / "kegg_escher_map.json").is_file()
    assert (output_dir / "merged_escher_map.json").is_file()
    assert (output_dir / "merged_kegg_escher_map.json").is_file()
    _assert_valid_saved_escher_map(output_dir / "merged_escher_map.json")
    _assert_valid_saved_escher_map(output_dir / "merged_kegg_escher_map.json")
    assert "merged_escher_map_json:" in completed.stdout
    assert "merged_kegg_escher_map_json:" in completed.stdout
