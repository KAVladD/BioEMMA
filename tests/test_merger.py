from copy import deepcopy

from bioemma.merger import EscherMerger
from bioemma.workflow import validate_escher_map


def _tiny_map(x_offset=0):
    return [
        {
            "map_name": "tiny",
            "map_id": "tiny",
            "map_description": "",
            "homepage": "https://escher.github.io",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#",
        },
        {
            "nodes": {
                "0": {"node_type": "metabolite", "x": x_offset, "y": 0},
                "1": {"node_type": "midmarker", "x": x_offset + 10, "y": 0},
            },
            "reactions": {
                "0": {
                    "bigg_id": "R_TEST",
                    "label_x": x_offset + 5,
                    "label_y": 0,
                    "segments": {
                        "0": {
                            "from_node_id": "0",
                            "to_node_id": "1",
                            "b1": None,
                            "b2": None,
                        }
                    },
                }
            },
            "text_labels": {},
            "canvas": {"x": 0, "y": 0, "width": 10, "height": 10},
        },
    ]


def test_merger_offsets_node_ids_and_keeps_segments_valid():
    first = _tiny_map()
    second = _tiny_map(100)

    merged = EscherMerger([deepcopy(first), deepcopy(second)]).merge_maps()

    assert len(merged[1]["nodes"]) == 4
    assert len(merged[1]["reactions"]) == 2
    validation = validate_escher_map(merged, strict_json_keys=True)
    assert validation["bad_segment_refs"] == []
    assert validation["duplicate_segment_ids"] == []
    assert validation["missing_description_keys"] == []
    assert validation["missing_model_keys"] == []

    second_reaction = merged[1]["reactions"]["1"]
    segment = next(iter(second_reaction["segments"].values()))
    assert next(iter(second_reaction["segments"].keys())) == "1"
    assert segment["from_node_id"] == "2"
    assert segment["to_node_id"] == "3"
    assert merged[1]["canvas"]["width"] > 0
    assert merged[1]["canvas"]["height"] > 0
