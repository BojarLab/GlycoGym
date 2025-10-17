import pandas as pd

from glycogym.glycogym import build_glycosylation, build_lgi, build_taxonomy, build_tissue


def test_linkage():
    df, mapping = build_glycosylation(top_k=1000)
    assert isinstance(df, pd.DataFrame)
    assert all(a == b for a, b in zip(df.columns, ["IUPAC", "SMILES", "label", "split"]))
    assert isinstance(mapping, dict)
    assert all(isinstance(k, str) and isinstance(v, int) for k, v in mapping.items())


def test_taxonomy():
    df = build_taxonomy("Kingdom", top_k=100)
    assert isinstance(df, pd.DataFrame)


def test_tissue():
    df = build_tissue(top_k=1000)
    assert isinstance(df, pd.DataFrame)


def test_lgi():
    df_r, df_cl, df_cg = build_lgi(top_k=100)
    assert isinstance(df_r, pd.DataFrame)
    assert isinstance(df_cl, pd.DataFrame)
    assert isinstance(df_cg, pd.DataFrame)
