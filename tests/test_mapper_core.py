from pathlib import Path

from bioemma.mapper_base import EscherMapper
from bioemma.workflow import build_outputs, validate_escher_map


ROOT = Path(__file__).resolve().parents[1]
KGML = ROOT / "tests" / "data" / "kgml" / "rn00010.xml"
MODEL = ROOT / "tests" / "data" / "models" / "e_coli_core.xml"


def test_build_kegg_map_keeps_first_reaction_index():
    metabolites = {
        "C00001": {
            "ids": {"KEGG": "C00001", "BIGG": "h2o", "SEED": "cpd00001"},
            "position": ("0", "0"),
        },
        "C00002": {
            "ids": {"KEGG": "C00002", "BIGG": "atp", "SEED": "cpd00002"},
            "position": ("100", "0"),
        },
    }
    reactions = {
        "first": {
            "ids": {"KEGG": "R00001", "BIGG": "FIRST", "SEED": "rxn00001"},
            "position": ("50", "0"),
            "substrates": {"main": ["C00001"], "side": []},
            "products": {"main": ["C00002"], "side": []},
            "reversibility": "irreversible",
        },
        "second": {
            "ids": {"KEGG": "R00002", "BIGG": "SECOND", "SEED": "rxn00002"},
            "position": ("150", "0"),
            "substrates": {"main": ["C00002"], "side": []},
            "products": {"main": ["C00001"], "side": []},
            "reversibility": "reversible",
        },
    }

    escher_map = EscherMapper(metabolites, reactions).build_kegg_map()

    assert len(escher_map[1]["reactions"]) == 2
    assert 0 in escher_map[1]["reactions"]
    assert escher_map[1]["reactions"][0]["name"] == "first"


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
