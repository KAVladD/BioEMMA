from pathlib import Path

from bioemma.maps import KeggMap


ROOT = Path(__file__).resolve().parents[1]
KGML = ROOT / "tests" / "data" / "kgml" / "rn00010.xml"


def _load_fixture_map() -> KeggMap:
    kegg_map = KeggMap()
    kegg_map.read_from_file(KGML.read_text(encoding="utf-8"))
    return kegg_map


def test_kegg_map_parses_compounds_reactions_and_details():
    kegg_map = _load_fixture_map()

    assert len(kegg_map.metabolites) == 31
    assert len(kegg_map.reactions) == 56

    metabolites = kegg_map.get_metabolites()
    reactions = kegg_map.get_reactions()

    assert metabolites["C00022"] == {
        "ids": {"KEGG": "C00022", "BIGG": "pyr", "SEED": "cpd00020"},
        "position": ("483", "868"),
    }
    assert reactions["R00200"]["ids"] == {
        "KEGG": "R00200",
        "BIGG": "PYK",
        "SEED": "rxn00148",
    }
    assert reactions["R00200"]["position"] == ("483", "771")
    assert reactions["R00200"]["substrates"] == {"main": ["C00074"], "side": []}
    assert reactions["R00200"]["products"] == {"main": ["C00022"], "side": []}
    assert reactions["R00200"]["reversibility"] == "irreversible"


def test_kegg_map_read_from_url_reads_response_without_temp_file(monkeypatch, tmp_path):
    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def read(self):
            return KGML.read_bytes()

    def fake_urlopen(url):
        assert url == "https://example.test/rn00010/kgml"
        return FakeResponse()

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("bioemma.maps.urlopen", fake_urlopen)

    kegg_map = KeggMap()
    kegg_map.read_from_url("https://example.test/rn00010/kgml")

    assert len(kegg_map.metabolites) == 31
    assert len(kegg_map.reactions) == 56
    assert not (tmp_path / "temp.xml").exists()
