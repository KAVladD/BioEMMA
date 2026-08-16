from pathlib import Path

import numpy as np

from bioemma.mapper_base import EscherMapper
from bioemma.metanetx_mapper import MetaNetXMapper
from bioemma.workflow import build_outputs, validate_escher_map


ROOT = Path(__file__).resolve().parents[1]
KGML = ROOT / "tests" / "data" / "kgml" / "rn00010.xml"
MODEL = ROOT / "tests" / "data" / "models" / "e_coli_core.xml"


def test_reaction_mapper_distinguishes_bigg_and_seed_ec_fallback_markers(tmp_path):
    mapping_path = tmp_path / "reaction_mapping.tsv"
    mapping_path.write_text(
        "\t".join(
            [
                "kegg",
                "mnx_id",
                "ec",
                "bigg",
                "seed",
                "metacyc",
                "rhea",
                "description",
                "ambiguous",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "R_SEED",
                "MNXR_SEED",
                "1.1.1.1",
                "",
                "rxnSEED",
                "",
                "",
                "",
                "seed_ec_fallback(participants>=0.50)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mapper = MetaNetXMapper(str(mapping_path))
    entry = mapper["R_SEED"]

    assert entry.is_seed_ec_fallback is True
    assert entry.is_bigg_ec_fallback is False
    assert entry.is_ec_fallback is False


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


def test_model_bigg_match_overrides_default_reaction_alias():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        def __init__(self, annotation):
            self.annotation = annotation
            self.metabolites = {}

    class FakeModel:
        def __init__(self, reaction):
            self.reactions = [reaction]
            self.metabolites = []

    bigg_model = FakeModel(FakeReaction({"bigg.reaction": "PGM"}))
    bigg_map = EscherMapper({}, reactions, scaling_factor=1).build_map(bigg_model)
    bigg_reaction = next(iter(bigg_map[1]["reactions"].values()))

    kegg_model = FakeModel(FakeReaction({"kegg.reaction": "R01518"}))
    kegg_map = EscherMapper({}, reactions, scaling_factor=1).build_map(kegg_model)
    kegg_reaction = next(iter(kegg_map[1]["reactions"].values()))

    assert bigg_reaction["bigg_id"] == "PGM"
    assert kegg_reaction["bigg_id"] == "PGAM_h"


def test_model_rhea_annotation_matches_reaction():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        def __init__(self, annotation):
            self.annotation = annotation
            self.metabolites = {}

    class FakeModel:
        reactions = [FakeReaction({"rhea": "15902"})]
        metabolites = []

    rhea_map = EscherMapper({}, reactions, scaling_factor=1).build_map(FakeModel())
    rhea_reaction = next(iter(rhea_map[1]["reactions"].values()))

    assert rhea_reaction["bigg_id"] == "PGAM_h"


def test_model_ec_annotation_matches_when_mapping_has_ec():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        def __init__(self, annotation):
            self.annotation = annotation
            self.metabolites = {}

    class FakeModel:
        reactions = [FakeReaction({"ec-code": "5.4.2.1"})]
        metabolites = []

    mapper = EscherMapper({}, reactions, scaling_factor=1)
    mapper.r_mapper["R01518"]._ec_all = ["5.4.2.1"]
    ec_map = mapper.build_map(FakeModel())
    ec_reaction = next(iter(ec_map[1]["reactions"].values()))

    assert ec_reaction["bigg_id"] == "PGAM_h"
    assert mapper.map_stats["model_matching"]["reaction_match_methods"] == {"ec": 1}
    assert mapper.map_stats["model_matching"]["use_fallback_matching"] is True


def test_model_ec_annotation_is_ignored_without_fallback_matching():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        annotation = {"ec-code": "5.4.2.1"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    mapper = EscherMapper({}, reactions, scaling_factor=1, use_fallback_matching=False)
    mapper.r_mapper["R01518"]._ec_all = ["5.4.2.1"]
    ec_map = mapper.build_map(FakeModel())

    assert ec_map[1]["reactions"] == {}
    assert mapper.map_stats["model_matching"]["matched_reactions"] == 0
    assert mapper.map_stats["model_matching"]["unmatched_reactions"] == 1
    assert mapper.map_stats["model_matching"]["reaction_match_methods"] == {}
    assert mapper.map_stats["model_matching"]["use_fallback_matching"] is False


def test_ec_fallback_bigg_mapping_is_ignored_without_fallback_matching():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        annotation = {"bigg.reaction": "PGM"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    strict_mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        use_fallback_matching=False,
    )
    strict_mapper.r_mapper["R01518"]._bigg_all = ["PGM"]
    strict_mapper.r_mapper["R01518"]._is_ec_fallback = True
    strict_map = strict_mapper.build_map(FakeModel())

    fallback_mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        use_fallback_matching=True,
        use_bigg_fallback_matching=True,
    )
    fallback_mapper.r_mapper["R01518"]._bigg_all = ["PGM"]
    fallback_mapper.r_mapper["R01518"]._is_ec_fallback = True
    fallback_map = fallback_mapper.build_map(FakeModel())

    assert strict_map[1]["reactions"] == {}
    fallback_reaction = next(iter(fallback_map[1]["reactions"].values()))
    assert fallback_reaction["bigg_id"] == "PGM"
    assert fallback_mapper.map_stats["model_matching"]["reaction_match_methods"] == {
        "bigg_ec_fallback": 1
    }


def test_bigg_ec_fallback_can_be_disabled_independently():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        annotation = {"bigg.reaction": "PGM"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        use_fallback_matching=True,
        use_bigg_fallback_matching=False,
        use_seed_fallback_matching=True,
    )
    mapper.r_mapper["R01518"]._bigg_all = ["PGM"]
    mapper.r_mapper["R01518"]._is_ec_fallback = True

    escher_map = mapper.build_map(FakeModel())

    assert escher_map[1]["reactions"] == {}
    assert mapper.map_stats["model_matching"]["reaction_match_methods"] == {}
    assert (
        mapper.map_stats["model_matching"]["use_bigg_fallback_matching"] is False
    )
    assert (
        mapper.map_stats["model_matching"]["use_seed_fallback_matching"] is True
    )


def test_seed_ec_fallback_can_be_disabled_independently():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "", "SEED": "rxnFALLBACK"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        id = "MODEL_SEED"
        annotation = {"seed.reaction": "rxnFALLBACK"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    disabled_mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        use_fallback_matching=True,
        use_seed_fallback_matching=False,
    )
    disabled_mapper.r_mapper["R01518"]._bigg_all = []
    disabled_mapper.r_mapper["R01518"]._seed_all = ["rxnFALLBACK"]
    disabled_mapper.r_mapper["R01518"]._is_seed_ec_fallback = True
    disabled_map = disabled_mapper.build_map(FakeModel())

    enabled_mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        use_fallback_matching=True,
        use_seed_fallback_matching=True,
    )
    enabled_mapper.r_mapper["R01518"]._bigg_all = []
    enabled_mapper.r_mapper["R01518"]._seed_all = ["rxnFALLBACK"]
    enabled_mapper.r_mapper["R01518"]._is_seed_ec_fallback = True
    enabled_map = enabled_mapper.build_map(FakeModel())

    assert disabled_map[1]["reactions"] == {}
    enabled_reaction = next(iter(enabled_map[1]["reactions"].values()))
    assert enabled_reaction["bigg_id"] == "MODEL_SEED"
    assert enabled_mapper.map_stats["model_matching"]["reaction_match_methods"] == {
        "seed_ec_fallback": 1
    }


def test_model_bigg_match_is_preferred_over_earlier_ec_match():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        def __init__(self, annotation):
            self.annotation = annotation
            self.metabolites = {}

    class FakeModel:
        reactions = [
            FakeReaction({"ec-code": "5.4.2.1"}),
            FakeReaction({"bigg.reaction": "PGM"}),
        ]
        metabolites = []

    mapper = EscherMapper({}, reactions, scaling_factor=1)
    mapper.r_mapper["R01518"]._ec_all = ["5.4.2.1"]
    mapper.r_mapper["R01518"]._bigg_all.append("PGM")

    escher_map = mapper.build_map(FakeModel())
    reaction = next(iter(escher_map[1]["reactions"].values()))

    assert reaction["bigg_id"] == "PGM"


def test_compartment_filter_selects_requested_model_reaction_compartment():
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
        "R00001": {
            "ids": {"KEGG": "R00001", "BIGG": "RXN_ONE", "SEED": "rxn00001"},
            "position": ("50", "0"),
            "substrates": {"main": ["C00001"], "side": []},
            "products": {"main": ["C00002"], "side": []},
            "reversibility": "irreversible",
        },
        "R00002": {
            "ids": {"KEGG": "R00002", "BIGG": "RXN_TWO", "SEED": "rxn00002"},
            "position": ("150", "0"),
            "substrates": {"main": ["C00001"], "side": []},
            "products": {"main": ["C00002"], "side": []},
            "reversibility": "irreversible",
        },
    }

    class FakeMetabolite:
        def __init__(self, metabolite_id, compartment, kegg_id, bigg_id):
            self.id = metabolite_id
            self.name = metabolite_id
            self.compartment = compartment
            self.annotation = {
                "kegg.compound": kegg_id,
                "bigg.metabolite": bigg_id,
            }

    class FakeReaction:
        def __init__(self, reaction_id, compartment):
            self.annotation = {"kegg.reaction": reaction_id}
            self.compartments = {compartment}
            h2o = FakeMetabolite(
                f"h2o_{compartment}",
                compartment,
                "C00001",
                "h2o",
            )
            atp = FakeMetabolite(
                f"atp_{compartment}",
                compartment,
                "C00002",
                "atp",
            )
            self.metabolites = {h2o: -1, atp: 1}

    class FakeModel:
        def __init__(self):
            self.reactions = [
                FakeReaction("R00001", "c"),
                FakeReaction("R00001", "m"),
                FakeReaction("R00002", "c"),
            ]
            self.metabolites = [
                metabolite
                for reaction in self.reactions
                for metabolite in reaction.metabolites
            ]

    mapper = EscherMapper(
        metabolites,
        reactions,
        scaling_factor=1,
        use_model_metabolite_ids=True,
        metabolite_id_compartments=True,
        compartment_filter=" m ",
    )

    escher_map = mapper.build_map(FakeModel())
    model = escher_map[1]
    reaction_names = {
        reaction["name"] for reaction in model["reactions"].values()
    }
    primary_nodes = [
        node
        for node in model["nodes"].values()
        if node.get("node_type") == "metabolite" and node.get("node_is_primary")
    ]

    assert reaction_names == {"R00001"}
    assert {node["bigg_id"] for node in primary_nodes} == {"h2o_m", "atp_m"}
    assert {node.get("compartment") for node in primary_nodes} == {"m"}
    assert mapper.map_stats["model_matching"]["compartment_filter"] == "m"
    assert mapper.map_stats["model_matching"]["matched_reactions"] == 1
    assert mapper.map_stats["model_matching"]["unmatched_reactions"] == 1


def test_non_bigg_model_match_uses_selected_model_reaction_id():
    reactions = {
        "R01518": {
            "ids": {"KEGG": "R01518", "BIGG": "PGAM_h", "SEED": "rxn01106"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        id = "PGM_c"
        annotation = {"ec-code": "5.4.2.1", "bigg.reaction": "PGM_c"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    mapper = EscherMapper({}, reactions, scaling_factor=1)
    mapper.r_mapper["R01518"]._ec_all = ["5.4.2.1"]
    mapper.r_mapper["R01518"]._bigg_all = ["PGAM_h"]

    escher_map = mapper.build_map(FakeModel())
    reaction = next(iter(escher_map[1]["reactions"].values()))

    assert reaction["bigg_id"] == "PGM_c"
    assert mapper.map_stats["model_matching"]["reaction_match_methods"] == {"ec": 1}


def test_duplicate_kegg_matches_to_same_model_reaction_keep_stronger_match():
    reactions = {
        "R01061": {
            "ids": {"KEGG": "R01061", "BIGG": "GAPD", "SEED": "rxn00781"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
        "R01063": {
            "ids": {"KEGG": "R01063", "BIGG": "GAPDH_nadp_hi", "SEED": "rxn00782"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        id = "GAPD"
        annotation = {"bigg.reaction": "GAPD", "ec-code": "1.2.1.12"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        use_bigg_fallback_matching=True,
    )
    mapper.r_mapper["R01061"]._bigg_all = ["GAPD"]
    mapper.r_mapper["R01063"]._bigg_all = ["GAPDH_nadp_hi"]
    mapper.r_mapper["R01063"]._ec_all = ["1.2.1.12"]

    escher_map = mapper.build_map(FakeModel())
    model_reactions = escher_map[1]["reactions"].values()

    assert {reaction["name"] for reaction in model_reactions} == {"R01061"}
    assert mapper.map_stats["model_matching"]["matched_reactions"] == 1
    assert mapper.map_stats["model_matching"]["unmatched_reactions"] == 1
    assert mapper.map_stats["model_matching"]["reaction_match_methods"] == {
        "bigg_ec_fallback": 1
    }


def test_distinct_model_reactions_at_same_kegg_position_are_kept():
    reactions = {
        "R01061": {
            "ids": {"KEGG": "R01061", "BIGG": "GAPD", "SEED": "rxn00781"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
        "R01063": {
            "ids": {"KEGG": "R01063", "BIGG": "GAPDy", "SEED": "rxn00782"},
            "position": ("50", "0"),
            "substrates": {"main": [], "side": []},
            "products": {"main": [], "side": []},
            "reversibility": "reversible",
        },
    }

    class FakeReaction:
        def __init__(self, reaction_id, ec):
            self.id = reaction_id
            self.annotation = {"bigg.reaction": reaction_id, "ec-code": ec}
            self.metabolites = {}

    class FakeModel:
        reactions = [
            FakeReaction("GAPD", "1.2.1.12"),
            FakeReaction("GAPDy", "1.2.1.13"),
        ]
        metabolites = []

    mapper = EscherMapper(
        {},
        reactions,
        scaling_factor=1,
        markers_dist=10,
        use_bigg_fallback_matching=True,
    )
    mapper.r_mapper["R01061"]._bigg_all = ["GAPD"]
    mapper.r_mapper["R01063"]._bigg_all = ["GAPDy"]

    escher_map = mapper.build_map(FakeModel())
    positions = {
        reaction["bigg_id"]: (float(reaction["label_x"]), float(reaction["label_y"]))
        for reaction in escher_map[1]["reactions"].values()
    }

    assert set(positions) == {"GAPD", "GAPDy"}


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


def test_remove_orphan_metabolites_removes_all_unreferenced_metabolite_nodes():
    mapper = EscherMapper({}, {})
    all_nodes = {
        1: {"node_type": "metabolite", "bigg_id": "primary", "node_is_primary": True},
        2: {"node_type": "metabolite", "bigg_id": "secondary", "node_is_primary": False},
        3: {"node_type": "metabolite", "bigg_id": "connected", "node_is_primary": True},
        4: {"node_type": "midmarker"},
    }
    r_desc = {
        "visible": {
            "segments": {
                0: {
                    "from_node_id": 3,
                    "to_node_id": 4,
                }
            }
        }
    }
    r2indx_dict = {"visible": 0}

    mapper._remove_orphan_metabolites(all_nodes, r_desc, r2indx_dict)

    assert all_nodes[1] is None
    assert all_nodes[2] is None
    assert all_nodes[3] is not None


def test_unmatched_metabolite_filter_keeps_nodes_referenced_by_visible_reactions():
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
        "R00001": {
            "ids": {"KEGG": "R00001", "BIGG": "RXN_ONE", "SEED": "rxn00001"},
            "position": ("50", "0"),
            "substrates": {"main": ["C00001"], "side": []},
            "products": {"main": ["C00002"], "side": []},
            "reversibility": "irreversible",
        },
    }

    class FakeReaction:
        annotation = {"kegg.reaction": "R00001"}
        metabolites = {}

    class FakeModel:
        reactions = [FakeReaction()]
        metabolites = []

    escher_map = EscherMapper(metabolites, reactions, scaling_factor=1).build_map(
        FakeModel()
    )
    validation = validate_escher_map(escher_map)
    primary_nodes = [
        node
        for node in escher_map[1]["nodes"].values()
        if node.get("node_type") == "metabolite" and node.get("node_is_primary")
    ]

    assert validation["bad_segment_refs"] == []
    assert {node["name"] for node in primary_nodes} == {"C00001", "C00002"}


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
    assert len(secondary_nodes) >= 26
    assert {"nad_c", "nadh_c", "adp_c"} <= {node["bigg_id"] for node in secondary_nodes}


def test_secondary_metabolites_can_use_selected_database_ids(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="SEED",
        scaling_factor=5,
        axis_epsilon=10,
        use_database_secondary_metabolite_ids=True,
    )

    nodes = result.escher_map[1]["nodes"]
    secondary_ids = {
        node["bigg_id"]
        for node in nodes.values()
        if node.get("node_type") == "metabolite" and not node.get("node_is_primary")
    }

    assert "cpd00003" in secondary_ids
    assert "nad_c" not in secondary_ids
    assert result.summary["metabolite_id_options"][
        "use_database_secondary_metabolite_ids"
    ]


def test_build_map_can_use_compartmental_model_metabolite_ids(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
        use_model_metabolite_ids=True,
        metabolite_id_compartments=True,
    )

    nodes = result.escher_map[1]["nodes"]
    primary_nodes = [
        node
        for node in nodes.values()
        if node.get("node_type") == "metabolite" and node.get("node_is_primary")
    ]
    secondary_nodes = [
        node
        for node in nodes.values()
        if node.get("node_type") == "metabolite" and not node.get("node_is_primary")
    ]

    glucose = next(node for node in primary_nodes if node["name"] == "C00031")
    assert glucose["bigg_id"] == "glc__D_e"
    assert glucose["compartment"] == "e"
    assert {"nad_c", "nadh_c", "adp_c"} <= {node["bigg_id"] for node in secondary_nodes}
    assert all("compartment" in node for node in secondary_nodes)


def test_build_map_can_use_model_metabolite_ids_without_compartments(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="BIGG",
        scaling_factor=5,
        axis_epsilon=10,
        use_model_metabolite_ids=True,
        metabolite_id_compartments=False,
    )

    nodes = result.escher_map[1]["nodes"]
    secondary_ids = {
        node["bigg_id"]
        for node in nodes.values()
        if node.get("node_type") == "metabolite" and not node.get("node_is_primary")
    }

    assert {"nad", "nadh", "adp"} <= secondary_ids
    assert {"nad_c", "nadh_c", "adp_c"}.isdisjoint(secondary_ids)
    assert all(
        "compartment" not in node
        for node in nodes.values()
        if node.get("node_type") == "metabolite"
    )


def test_seed_model_ids_keep_database_ids_and_store_compartments(monkeypatch, tmp_path):
    cache_root = tmp_path / "cobra-cache"
    monkeypatch.setenv("BIOEMMA_COBRA_CACHE_DIR", str(cache_root))

    result = build_outputs(
        model=MODEL,
        kgml=KGML,
        database="SEED",
        scaling_factor=5,
        axis_epsilon=10,
        use_model_metabolite_ids=True,
        metabolite_id_compartments=True,
    )

    nodes = result.escher_map[1]["nodes"]
    secondary_nodes = [
        node
        for node in nodes.values()
        if node.get("node_type") == "metabolite" and not node.get("node_is_primary")
    ]

    assert "cpd00003" in {node["bigg_id"] for node in secondary_nodes}
    assert "nad_c" not in {node["bigg_id"] for node in secondary_nodes}
    assert all(not node["bigg_id"].endswith("_c") for node in secondary_nodes)
    assert all(node.get("compartment") == "c" for node in secondary_nodes)


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
