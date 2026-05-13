from bioemma import EscherMapper, KeggMap, MetaNetXMapper


def test_public_imports():
    assert KeggMap is not None
    assert EscherMapper is not None
    assert MetaNetXMapper is not None
