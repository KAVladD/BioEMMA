from pathlib import Path

from bioemma.workflow import build_outputs, validate_escher_map


ROOT = Path(__file__).resolve().parents[1]
KGML = ROOT / "tests" / "data" / "kgml" / "rn00010.xml"
MODEL = ROOT / "tests" / "data" / "models" / "e_coli_core.xml"


def test_build_map_adds_secondary_metabolites_and_valid_segments(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
    )

    validation = validate_escher_map(result.escher_map)
    assert validation["bad_segment_refs"] == []

    nodes = result.escher_map[1]["nodes"]
    secondary_nodes = [
        node
        for node in nodes.values()
        if node.get("node_type") == "metabolite" and not node.get("node_is_primary")
    ]
    assert len(secondary_nodes) == 26
    assert {"nad_c", "nadh_c", "adp_c"} <= {node["bigg_id"] for node in secondary_nodes}


def test_build_map_canvas_encloses_generated_nodes(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
    )

    model = result.escher_map[1]
    canvas = model["canvas"]
    xs = [float(node["x"]) for node in model["nodes"].values()]
    ys = [float(node["y"]) for node in model["nodes"].values()]

    assert canvas["width"] > 0
    assert canvas["height"] > 0
    assert min(xs) >= canvas["x"]
    assert min(ys) >= canvas["y"]
    assert max(xs) <= canvas["x"] + canvas["width"]
    assert max(ys) <= canvas["y"] + canvas["height"]
