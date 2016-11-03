//
// Created by M. Engler on 26/10/16.
//

#include <chrono>
#include "product.h"

using namespace mogli;

void Product::dfs(Molecule &mol, const Node &v, int depth, NodeToBoolMap &visited, NodeToBoolMap &filter) {
  visited[v] = true;
  filter[v] = true;
  if (depth < _shell) {
    for (IncEdgeIt e = mol.get_inc_edge_iter(v); e != lemon::INVALID; ++e) {
      Node w = mol.get_opposite_node(v, e);
      if (!visited[w]) {
        dfs(mol, w, depth+1, visited, filter);
      }
    }
  }
}

void Product::generate_subgraph_canonization(Molecule &mol, const Node &v, NodeToCanonizationMap &map) {
  NodeToBoolMap filter(mol.get_graph(), false);
  NodeToBoolMap visited(mol.get_graph(), false);
  dfs(mol, v, 0, visited, filter);
  map[v] = Canonization(mol, filter);
};

void Product::generate0() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  // generate nodes
  for (NodeIt u = _mol1.get_node_iter(); u != lemon::INVALID; ++u) {
    for (NodeIt v = _mol2.get_node_iter(); v != lemon::INVALID; ++v) {
      if (_mol1.get_color(u) == _mol2.get_color(v)) {
        Node uv = _g.addNode();
        _mol1ToG[u] = uv;
        _mol2ToG[v] = uv;
        _gToMol1[uv] = u;
        _gToMol2[uv] = v;
        ++_numNodes;
      }
    }
  }

  // generate edges
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
    Node u1 = _gToMol1[u1v1];
    Node v1 = _gToMol2[u1v1];
    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2) {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _gToMol1[u2v2];
      Node v2 = _gToMol2[u2v2];

      assert(_mol1.get_color(u1) == _mol2.get_color(v1));

      if (u1 != u2 && v1 != v2) {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 == v1v2) {
          _connectivityEdge[_g.addEdge(u1v1, u2v2)] = u1u2;
          ++_numEdges;
        }
      }
    }
  }

}

void Product::generate() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  NodeToCanonizationMap subgraph_canons1(g1);
  NodeToCanonizationMap subgraph_canons2(g2);

  for (NodeIt v(g1); v != lemon::INVALID; ++v) {
    generate_subgraph_canonization(_mol1, v, subgraph_canons1);
  }
  for (NodeIt v(g2); v != lemon::INVALID; ++v) {
    generate_subgraph_canonization(_mol2, v, subgraph_canons2);
  }

  NodeMatrix matrix;
  matrix.reserve(_mol1.get_atom_count());

  // generate nodes
  for (NodeIt u = _mol1.get_node_iter(); u != lemon::INVALID; ++u) {
    NodeVector row;
    for (NodeIt v = _mol2.get_node_iter(); v != lemon::INVALID; ++v) {
      if (_mol1.get_color(u) == _mol2.get_color(v) &&
          are_isomorph(subgraph_canons1[u],subgraph_canons2[v])) {
        Node uv = _g.addNode();
//        std::cout << _g.id(uv) << " " << _mol1.get_label(u) << " " << _mol1.get_chem_element(u) << " " << _mol1.get_color(u) << " "
//            << " " << _mol2.get_label(v) << " " << _mol2.get_chem_element(v) << " " << _mol2.get_color(v) << std::endl;
        _mol1ToG[u] = uv;
        _mol2ToG[v] = uv;
        _gToMol1[uv] = u;
        _gToMol2[uv] = v;
        ++_numNodes;
        row.push_back(uv);
      }
    }
    matrix.push_back(row);
  }

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  for (NodeMatrix::const_iterator it = matrix.begin(), end = matrix.end(); it != end; ++it) {
    for (NodeVector::const_iterator it2 = (*it).begin(), end2 = (*it).end(); it2 != end2; ++it2) {
      Node u1v1 = (*it2);
      Node u1 = _gToMol1[u1v1];
      Node v1 = _gToMol2[u1v1];

      for (NodeMatrix::const_iterator it3 = it+1; it3 != end; ++it3) {
        for (NodeVector::const_iterator it4 = (*it3).begin(), end4 = (*it3).end(); it4 != end4; ++it4) {
          Node u2v2 = (*it4);
          Node u2 = _gToMol1[u2v2];
          Node v2 = _gToMol2[u2v2];

          if (u1 == u2 || v1 == v2)
            continue;

          bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
          bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

          if (u1u2 == v1v2) {
            _connectivityEdge[_g.addEdge(u1v1, u2v2)] = u1u2;
            ++_numEdges;
          }
        }
      }

    }
  }

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  std::cout << "new: " << duration << std::endl;

}

void Product::determineDegrees(const Graph& g, NodeToIntMap& deg) {
  for (NodeIt v(g); v != lemon::INVALID; ++v) {
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e) {
      ++deg[v];
    }
  }
}

void Product::generate(const Molecule& mol,
                                  const NodeToIntMap& deg,
                                  IntSetNodeMap& intSet,
                                  IntSetNodeMap& degSet) {
  const Graph& g = mol.get_graph();
  NodeToBoolMap visited(g, false);
  for (NodeIt v(g); v != lemon::INVALID; ++v)
  {
    lemon::mapFill(g, visited, false);
    IntSet& s = intSet[v];
    IntSet& ds = degSet[v];
    dfs(deg, v, 0, mol, visited, s, ds);
  }
}

void Product::dfs(const NodeToIntMap& deg,
                     const Node v, const int depth,
                     const Molecule& mol,
                     NodeToBoolMap& visited,
                     IntSet& s, IntSet& ds) {
  const Graph& g = mol.get_graph();
  visited[v] = true;
  s.insert(mol.get_color(v));
  ds.insert(deg[v]);

  if (depth < _shell)
  {
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e)
    {
      Node w = g.oppositeNode(v, e);
      if (!visited[w])
      {
        dfs(deg, w, depth + 1, mol, visited, s, ds);
      }
    }
  }
}

void Product::generateDeg1NeighborSet(const Graph& g,
                                           const NodeToIntMap& deg,
                                           NodeVectorMap& deg1NeighborMap) {
  for (NodeIt u(g); u != lemon::INVALID; ++u) {
    if (deg[u] > 1) {
      for (IncEdgeIt e(g, u); e != lemon::INVALID; ++e) {
        Node v = g.oppositeNode(u, e);
        if (deg[v] == 1) {
          deg1NeighborMap[u].push_back(v);
        }
      }
    }
  }
}

void Product::generate_old() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  NodeToIntSetMap set1(g1);
  NodeToIntSetMap degSet1(g1);
  NodeToIntSetMap set2(g2);
  NodeToIntSetMap degSet2(g2);

  NodeToIntMap deg1(g1, 0);
  determineDegrees(g1, deg1);
  NodeToIntMap deg2(g2, 0);
  determineDegrees(g2, deg2);

  // determine degrees
  generate(_mol1, deg1, set1, degSet1);
  generate(_mol2, deg2, set2, degSet2);

  // generate nodes
  for (NodeIt u(g1); u != lemon::INVALID; ++u) {
    for (NodeIt v(g2); v != lemon::INVALID; ++v) {
      if (_mol1.get_color(u) == _mol2.get_color(v) && set1[u] == set2[v] && degSet1[u] == degSet2[v]) {
        assert(deg1[u] == deg2[v]);
        // don't add product nodes for a pair of deg-1 nodes unless _shell == 0
//        if (_shell == 0 || deg1[u] > 1) {
          Node uv = _g.addNode();
          _mol1ToG[u] = uv;
          _mol2ToG[v] = uv;
          _gToMol1[uv] = u;
          _gToMol2[uv] = v;
          ++_numNodes;
//        }
      }
    }
  }

  // generate deg1 sets
  generateDeg1NeighborSet(g1, deg1, _g1ToDeg1Neighbors);
  generateDeg1NeighborSet(g2, deg2, _g2ToDeg1Neighbors);

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

  // generate edges
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
    Node u1 = _gToMol1[u1v1];
    Node v1 = _gToMol2[u1v1];
    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2) {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _gToMol1[u2v2];
      Node v2 = _gToMol2[u2v2];

      assert(_mol1.get_color(u1) == _mol2.get_color(v1));

      if (u1 != u2 && v1 != v2) {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 == v1v2) {
          _connectivityEdge[_g.addEdge(u1v1, u2v2)] = u1u2;
          ++_numEdges;
        }
      }
    }
  }

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

  std::cout << "old: " << duration << std::endl;
}
