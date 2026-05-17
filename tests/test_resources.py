from bioemma._resources import resource_path
from bioemma.metanetx_mapper import MetaNetXMapper


def test_bundled_mapping_resources_exist():
    metabolite_mapping = resource_path("metabolite_mapping.tsv")
    reaction_mapping = resource_path("reaction_mapping.tsv")

    assert metabolite_mapping.is_file()
    assert reaction_mapping.is_file()


def test_bundled_metabolite_mapping_can_be_loaded():
    mapper = MetaNetXMapper(resource_path("metabolite_mapping.tsv"))

    assert len(mapper) > 0
    assert "C00002" in mapper
    assert mapper["C00002"].bigg


def test_bundled_reaction_mapping_can_be_loaded():
    mapper = MetaNetXMapper(resource_path("reaction_mapping.tsv"))

    assert len(mapper) > 0
    assert "R00200" in mapper
    assert mapper["R00200"].bigg == "PYK"
    assert mapper["R00200"].seed == "rxn00148"
