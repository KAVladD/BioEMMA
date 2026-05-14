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


def test_build_outputs_returns_and_saves_core_artifacts(monkeypatch, tmp_path):
    _set_cobra_cache(monkeypatch, tmp_path)

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        output_dir=tmp_path / "out",
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
    )

    assert result.escher_map[0]["schema"] == "https://escher.github.io/escher/jsonschema/1-0-0#"
    assert result.kegg_reconstruction["counts"]["metabolites"] > 0
    assert result.kegg_reconstruction["counts"]["reactions"] > 0

    validation = validate_escher_map(result.escher_map)
    assert validation["nodes"] > 0
    assert validation["reactions"] > 0
    assert validation["segments"] > 0
    assert validation["bad_segment_refs"] == []

    assert result.paths["map_json"].is_file()
    assert result.paths["kegg_reconstruction_json"].is_file()
    assert result.paths["summary_json"].is_file()

    saved_map = json.loads(result.paths["map_json"].read_text(encoding="utf-8"))
    saved_reconstruction = json.loads(
        result.paths["kegg_reconstruction_json"].read_text(encoding="utf-8")
    )
    assert validate_escher_map(saved_map)["bad_segment_refs"] == []
    assert saved_reconstruction["counts"] == result.kegg_reconstruction["counts"]


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
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )

    output_dir = tmp_path / "cli-out" / "rn00010"
    assert (output_dir / "map.json").is_file()
    assert (output_dir / "kegg_reconstruction.json").is_file()
    assert (output_dir / "summary.json").is_file()
    assert "map_json:" in completed.stdout


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
    assert (tmp_path / "batch-out" / "rn00010" / "map.json").is_file()
    assert (tmp_path / "batch-out" / "rn00020" / "map.json").is_file()
    assert result.paths["merged_map_json"].is_file()

    validation = validate_escher_map(result.merged_map)
    assert validation["nodes"] > 0
    assert validation["reactions"] > 0
    assert validation["bad_segment_refs"] == []


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
        ],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )

    output_dir = tmp_path / "cli-batch-out"
    assert (output_dir / "rn00010" / "map.json").is_file()
    assert (output_dir / "rn00020" / "map.json").is_file()
    assert (output_dir / "merged_map.json").is_file()
    assert "merged_map_json:" in completed.stdout
