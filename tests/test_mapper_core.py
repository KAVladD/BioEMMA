from pathlib import Path

import numpy as np

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


def test_mapper_accepts_layout_visualization_options():
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
    }

    escher_map = EscherMapper(
        metabolites,
        reactions,
        scaling_factor=1,
        metabolite_label_shift=(3, 7),
        reaction_label_shift=(11, 13),
        canvas_margin_x=250,
        canvas_margin_y=125,
        axis_offset=15,
    ).build_kegg_map()

    model = escher_map[1]
    metabolite = next(
        node for node in model["nodes"].values() if node.get("bigg_id") == "h2o"
    )
    midmarker = next(
        node for node in model["nodes"].values() if node.get("node_type") == "midmarker"
    )
    reaction = model["reactions"][0]

    assert float(metabolite["label_x"]) - float(metabolite["x"]) == 3
    assert float(metabolite["label_y"]) - float(metabolite["y"]) == 7
    assert float(reaction["label_x"]) - float(midmarker["x"]) == 11
    assert float(reaction["label_y"]) - float(midmarker["y"]) == 13
    assert model["canvas"]["width"] == 350
    assert model["canvas"]["height"] == 125


def test_multimarker_distance_options_affect_aligned_reactions():
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
    }

    mapper = EscherMapper(
        metabolites,
        reactions,
        scaling_factor=1,
        multimarker_distance_fraction=0.5,
    )
    model = mapper.build_kegg_map()[1]
    midmarker = next(
        node for node in model["nodes"].values() if node.get("node_type") == "midmarker"
    )
    multimarkers = [
        node
        for node in model["nodes"].values()
        if node.get("node_type") == "multimarker"
    ]

    assert sorted(abs(node["x"] - midmarker["x"]) for node in multimarkers) == [25, 25]

    constant_mapper = EscherMapper(
        metabolites,
        reactions,
        scaling_factor=1,
        use_constant_multimarker_distance=True,
        constant_multimarker_distance=7,
    )
    constant_model = constant_mapper.build_kegg_map()[1]
    constant_midmarker = next(
        node
        for node in constant_model["nodes"].values()
        if node.get("node_type") == "midmarker"
    )
    constant_multimarkers = [
        node
        for node in constant_model["nodes"].values()
        if node.get("node_type") == "multimarker"
    ]

    assert sorted(
        abs(node["x"] - constant_midmarker["x"]) for node in constant_multimarkers
    ) == [7, 7]


def test_secondary_metabolite_spacing_options_affect_positions():
    default_mapper = EscherMapper({}, {}, markers_dist=5)
    inherited = default_mapper._calc_secondary_position(
        anchor_pos=np.array([10.0, 10.0]),
        reaction_dir=np.array([1.0, 0.0]),
        perp=np.array([0.0, 1.0]),
        index=1,
        total=3,
        side=1,
    )
    assert inherited == [20.0, 17.5]

    mapper = EscherMapper(
        {},
        {},
        secondary_metabolite_distance=12,
        secondary_metabolite_spacing=8,
    )

    single = mapper._calc_secondary_position(
        anchor_pos=np.array([10.0, 10.0]),
        reaction_dir=np.array([1.0, 0.0]),
        perp=np.array([0.0, 1.0]),
        index=0,
        total=1,
        side=1,
    )
    first = mapper._calc_secondary_position(
        anchor_pos=np.array([10.0, 10.0]),
        reaction_dir=np.array([1.0, 0.0]),
        perp=np.array([0.0, 1.0]),
        index=0,
        total=3,
        side=1,
    )
    last = mapper._calc_secondary_position(
        anchor_pos=np.array([10.0, 10.0]),
        reaction_dir=np.array([1.0, 0.0]),
        perp=np.array([0.0, 1.0]),
        index=2,
        total=3,
        side=1,
    )

    assert single == [22.0, 14.0]
    assert first == [22.0, 2.0]
    assert last == [22.0, 18.0]
    assert mapper._calc_secondary_lateral_offset(0, 1) != 0
    assert mapper._calc_secondary_lateral_offset(1, 3) != 0


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
