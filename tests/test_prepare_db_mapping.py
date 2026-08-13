from scripts.prepare_db_mapping import (
    build_ec_fallback_offline,
    build_seed_ec_fallback_offline,
    parse_reaction_participants,
    participant_overlap,
)


def test_parse_reaction_participants_extracts_unique_mnx_ids():
    equation = "1 MNXM1@MNXD1 + 2 WATER@MNXD2 = 1 MNXM2@MNXD1 + 1 MNXM1@MNXD3"

    assert parse_reaction_participants(equation) == ["MNXM1", "MNXM2", "WATER"]


def test_participant_overlap_uses_source_denominator():
    source = {"A", "B", "C", "D"}
    candidate = {"A", "B", "X"}

    assert participant_overlap(source, candidate) == 0.5
    assert participant_overlap(set(), candidate) == 0.0


def test_ec_fallback_filters_candidates_by_participant_overlap():
    mnx_mapping = {
        "source": {
            "kegg": ["R_SOURCE"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "B", "C", "D"],
        },
        "good": {
            "bigg": ["GOOD"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "B", "C", "X"],
        },
        "lower_overlap": {
            "bigg": ["LOWER"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "B", "X"],
        },
        "bad": {
            "bigg": ["BAD"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "Y"],
        },
    }
    rows = [
        {
            "kegg": "R_SOURCE",
            "mnx_id": "source",
            "ec": "1.1.1.1",
            "bigg": "",
            "ambiguous": "",
        }
    ]

    build_ec_fallback_offline(mnx_mapping, rows)

    assert rows[0]["bigg"] == "GOOD"
    assert rows[0]["ambiguous"] == "ec_fallback(participants>=0.50)"


def test_seed_ec_fallback_filters_candidates_by_participant_overlap():
    mnx_mapping = {
        "source": {
            "kegg": ["R_SOURCE"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "B", "C", "D"],
        },
        "good": {
            "seed": ["rxnGOOD"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "B", "C", "X"],
        },
        "lower_overlap": {
            "seed": ["rxnLOWER"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "B", "X"],
        },
        "bad": {
            "seed": ["rxnBAD"],
            "ec": ["1.1.1.1"],
            "participants": ["A", "Y"],
        },
    }
    rows = [
        {
            "kegg": "R_SOURCE",
            "mnx_id": "source",
            "ec": "1.1.1.1",
            "seed": "",
            "ambiguous": "",
        }
    ]

    build_seed_ec_fallback_offline(mnx_mapping, rows)

    assert rows[0]["seed"] == "rxnGOOD"
    assert rows[0]["ambiguous"] == "seed_ec_fallback(participants>=0.50)"
