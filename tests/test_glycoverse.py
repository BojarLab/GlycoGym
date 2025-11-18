from typing import Generator

import networkx
import pytest
from glycogym.glycoverse import get_all_iupacs, generate_minus_one_subgraphs, explore_subgraph_topology, explore_glycoverse
from glycowork.motif.graph import glycan_to_nxGraph, graph_to_string, get_possible_topologies

def test_get_all_iupacs():
    graph = glycan_to_nxGraph("Man(a1-4)[Glc(a1-3)]Man")
    iupacs = get_all_iupacs(graph)
    assert isinstance(iupacs, list)
    assert all(isinstance(i, str) for i in iupacs)
    assert "Man(a1-4)[Glc(a1-3)]Man" in iupacs
    assert "Glc(a1-3)[Man(a1-4)]Man" in iupacs


def test_generate_minus_one_subgraphs():
    graph = glycan_to_nxGraph("Man(a1-4)[Glc(a1-3)]Man")
    subgen = generate_minus_one_subgraphs(graph)
    subgraphs = list(subgen)
    iupacs = [graph_to_string(g) for g in subgraphs]

    assert isinstance(subgen, Generator)
    assert all(isinstance(g, type(graph)) for g in subgraphs)
    assert len(subgraphs) == 2  # Two leaf nodes to remove
    assert "Man(a1-4)Man" in iupacs
    assert "Glc(a1-3)Man" in iupacs


def test_explore_subgraph_topology():
    graph = glycan_to_nxGraph("Man(a1-4)[Glc(a1-3)]Man")
    unique_iupacs = explore_subgraph_topology(graph)

    assert isinstance(unique_iupacs, set)
    assert all(isinstance(i, str) for i in unique_iupacs)
    assert "Man(a1-4)Man" in unique_iupacs
    assert "Glc(a1-3)Man" in unique_iupacs
    assert "Man" in unique_iupacs


def test_get_possible_topologies():
    """Copied from gycowork to ensure their function works as a proxy to test get_all_topologies_ext."""
    glycan = "{Neu5Ac(a2-?)}Gal(b1-3)GalNAc"
    topologies = get_possible_topologies(glycan)
    topologies = get_possible_topologies(glycan, return_graphs=True)
    assert len(topologies) > 0
    assert all(isinstance(g, networkx.Graph) for g in topologies)
    assert get_possible_topologies("{6S}Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-6)]GalNAc")[0] == "Neu5Ac(a2-3)Gal6S(b1-3)[Neu5Ac(a2-6)]GalNAc"
    assert get_possible_topologies("{OS}Gal(b1-3)GalNAc") == ['GalOS(b1-3)GalNAc', 'Gal(b1-3)GalNAcOS']
    with pytest.raises(ValueError, match="This glycan already has a defined topology; please don't use this function."):
        get_possible_topologies("Gal(b1-4)GlcNAc")


def test_explore_glycoverse():
    iupacs = explore_glycoverse(["Man(a1-4)[Glc(a1-3)]Man", "Gal(b1-4)GlcNAc", "Neu5Ac(a2-6)Gal(b1-3)GlcNAc"])

    assert isinstance(iupacs, set)
    assert all(isinstance(i, str) for i in iupacs)
    assert len(iupacs) == 8
    assert "Man(a1-4)[Glc(a1-3)]Man" in iupacs
    assert "Man(a1-4)Man" in iupacs
    assert "Glc(a1-3)Man" in iupacs
    assert "Man" in iupacs
    assert "Gal(b1-4)GlcNAc" in iupacs
    assert "GlcNAc" in iupacs
    assert "Neu5Ac(a2-6)Gal(b1-3)GlcNAc" in iupacs
    assert "Gal(b1-3)GlcNAc" in iupacs
