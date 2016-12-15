//
// Created by M. Engler on 26/10/16.
//

#include "product.h"

using namespace mogli;

void Product::determine_degrees(const Graph& g, IntToNodeMap& deg_to_node, NodeToIntMap& deg) {
  for (NodeIt v(g); v != lemon::INVALID; ++v) {
    int degree = 0;
    for (IncEdgeIt e(g, v); e != lemon::INVALID; ++e) {
      ++degree;
    }
    deg_to_node.insert(std::make_pair(degree, v));
    deg[v] = degree;
  }
}

void Product::bfs(const Molecule &mol, const Node &v, NodeToBoolMap &filter) {
  NodeToBoolMap visited(mol.get_graph(), false);
  NodeToIntMap depth(mol.get_graph(), 0);
  NodeDeque queue;

  queue.push_back(v);
  visited[v] = true;
  while (queue.size() > 0) {
    Node &current = queue.front();
    filter[current] = true;

    if (depth[current] < _shell) {
      for (IncEdgeIt e = mol.get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
        Node w = mol.get_opposite_node(current, e);
        if (!visited[w]) {
          visited[w] = true;
          depth[w] = depth[current]+1;
          queue.push_back(w);
        }
      }
    }
    queue.pop_front();
  }
}

void Product::bfs_neighbors(const Molecule &mol, const Node &v, BitSet &neighbors, int& size) {
  NodeToBoolMap visited(mol.get_graph(), false);
  NodeToIntMap depth(mol.get_graph(), 0);
  NodeDeque queue;

  queue.push_back(v);
  visited[v] = true;
  while (queue.size() > 0) {
    Node &current = queue.front();
    neighbors[mol.get_id(current)] = true;
    ++size;

    if (depth[current] < _shell) {
      for (IncEdgeIt e = mol.get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
        Node w = mol.get_opposite_node(current, e);
        if (!visited[w]) {
          visited[w] = true;
          depth[w] = depth[current]+1;
          queue.push_back(w);
        }
      }
    }
    queue.pop_front();
  }
}

void Product::bfs_subgraph(const Molecule &mol, const Node &product_node, const Node &root_node, const BitSet &root_neighbors,
                           const NodeToBitSetMap &neighborhoods, const NodeVector &order1, const NodeVector &order2,
                           ShortToNodeVectorPairMap &current_reductions) {
  NodeToBoolMap visited(mol.get_graph(), false);
  NodeToIntMap depth(mol.get_graph(), 0);
  NodeDeque queue;

  queue.push_back(root_node);
  visited[root_node] = true;
  while (queue.size() > 0) {
    Node &current = queue.front();

    const BitSet &neighbors_current = neighborhoods[current];
    if (neighbors_current.is_subset_of(root_neighbors)) {
      if (current != root_node) {
        int index = static_cast<int>(std::find(order1.begin(), order1.end(), current) - order1.begin());
        Node w = order2[index];
        unsigned short color = _mol1.get_color(current);
        if (current_reductions.find(color) == current_reductions.end()) {
          current_reductions[color] = std::make_pair(NodeVector(), NodeVector());
        }
        current_reductions[color].first.push_back(current);
        current_reductions[color].second.push_back(w);
        _reductions[product_node].push_back(std::make_pair(current, w));
      }

      if (depth[current] < _shell) {
        for (IncEdgeIt e = mol.get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
          Node x = mol.get_opposite_node(current, e);
          if (!visited[x]) {
            visited[x] = true;
            depth[x] = depth[current] + 1;
            queue.push_back(x);
          }
        }
      }
    }
    queue.pop_front();
  }
}

void Product::generate_subgraph_canonization(const Molecule &mol, const Node &v, NodeToCanonizationMap &map) {
  NodeToBoolMap filter(mol.get_graph(), false);
  bfs(mol, v, filter);
  map[v] = Canonization(mol, filter, v);
}

void Product::generate_subgraph(const Molecule &mol, const Node &v, NodeToBitSetMap &neighborhoods, IntToNodeMap &sizes) {
  int size = 0;
  bfs_neighbors(mol, v, neighborhoods[v], size);
  sizes.insert(std::make_pair(size, v));
}

void Product::generate_nodes() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  NodeToCanonizationMap subgraph_canons1(g1);
  NodeToCanonizationMap subgraph_canons2(g2);

  for (NodeIt v(g1); v != lemon::INVALID; ++v) {
    generate_subgraph_canonization(_mol1, v, subgraph_canons1);
  }
  for (NodeIt v(g2); v != lemon::INVALID; ++v) {
    generate_subgraph_canonization(_mol2, v, subgraph_canons2);
  }

  // generate product nodes
  for (NodeIt u = _mol1.get_node_iter(); u != lemon::INVALID; ++u) {
    for (NodeIt v = _mol2.get_node_iter(); v != lemon::INVALID; ++v) {
      if (_mol1.get_color(u) == _mol2.get_color(v) &&
          are_isomorphic(subgraph_canons1[u], subgraph_canons2[v])) {
        const Node& uv = add_node(u, v);
      }
    }
  }
}

void Product::generate_nodes_deg1() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  IntToNodeMap deg_to_node1;
  IntToNodeMap deg_to_node2;

  NodeToIntMap deg1(g1);
  NodeToIntMap deg2(g2);

  determine_degrees(g1, deg_to_node1, deg1);
  determine_degrees(g2, deg_to_node2, deg2);

  NodeToCanonizationMap subgraph_canons(g2);
  NodeToBoolMap has_canon(g2, false);

  NodeToBitSetMap reduced_nodes(g1, BitSet(_mol2.get_atom_count()));

  // iterate nodes u with decreasing degree
  for (IntToNodeMap::reverse_iterator it = deg_to_node1.rbegin(), end = deg_to_node1.rend(); it != end; ++it) {
    Node u = it->second;
    NodeToBoolMap filter(g1, false);
    bfs(_mol1, u, filter);
    Canonization canon1 = Canonization(_mol1, filter, u);
    // iterate nodes v with decreasing degree
    for (IntToNodeMap::reverse_iterator it2 = deg_to_node2.rbegin(), end2 = deg_to_node2.rend(); it2 != end2; ++it2) {
      Node v = it2->second;
      if (reduced_nodes[u][_mol2.get_id(v)]) {
        continue;
      }
      if (_mol1.get_color(u) == _mol2.get_color(v)) {
        if (!has_canon[v]) {
          generate_subgraph_canonization(_mol2, v, subgraph_canons);
          has_canon[v] = true;
        }
        Canonization& canon2 = subgraph_canons[v];
        if (are_isomorphic(canon1, canon2)) {
          // generate product node uv
          const Node& uv = add_node(u, v);
          // apply degree-1 rule to neighbors of u and v
          const NodeVector& order1 = canon1.get_node_order();
          const NodeVector& order2 = canon2.get_node_order();

          ShortToNodeVectorPairMap current_reductions;
          for (IncEdgeIt e(g1, u); e != lemon::INVALID; ++e) {
            Node w = g1.oppositeNode(u, e);
            if (deg1[w] == 1) {
              int index = static_cast<int>(std::find(order1.begin(), order1.end(), w) - order1.begin());
              Node x = order2[index];
              unsigned short color = _mol1.get_color(w);
              if (current_reductions.find(color) == current_reductions.end()) {
                current_reductions[color] = std::make_pair(NodeVector(), NodeVector());
              }
              current_reductions[color].first.push_back(w);
              current_reductions[color].second.push_back(x);
              _reductions[uv].push_back(std::make_pair(w, x));
            }
          }
          for (ShortToNodeVectorPairMap::const_iterator it3 = current_reductions.begin(), end3 = current_reductions.end(); it3 != end3; ++it3) {
            for (int i = 0; i < it3->second.first.size(); ++i) {
              for (int j = 0; j < it3->second.second.size(); ++j) {
                reduced_nodes[it3->second.first[i]][_mol2.get_id(it3->second.second[j])] = true;
              }
            }
          }
        }
      }
    }
  }

}

void Product::generate_nodes_sub() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  IntToNodeMap neighborhood_sizes1;
  IntToNodeMap neighborhood_sizes2;

  NodeToBitSetMap neighborhoods1(g1, BitSet(_mol1.get_atom_count()));
  NodeToBitSetMap neighborhoods2(g2, BitSet(_mol2.get_atom_count()));

  for (NodeIt v(g1); v != lemon::INVALID; ++v) {
    generate_subgraph(_mol1, v, neighborhoods1, neighborhood_sizes1);
  }

  for (NodeIt v(g2); v != lemon::INVALID; ++v) {
    generate_subgraph(_mol2, v, neighborhoods2, neighborhood_sizes2);
  }

  NodeToCanonizationMap subgraph_canons(g2);
  NodeToBoolMap has_canon(g2, false);

  NodeToBitSetMap reduced_nodes(g1, BitSet(_mol2.get_atom_count()));

  // iterate nodes u with decreasing neighborhood size
  for (IntToNodeMap::reverse_iterator it = neighborhood_sizes1.rbegin(), end = neighborhood_sizes1.rend(); it != end; ++it) {
    Node u = it->second;
    NodeToBoolMap filter(g1, false);
    bfs(_mol1, u, filter);
    Canonization canon1 = Canonization(_mol1, filter, u);
    // iterate nodes v with decreasing neighborhood size
    for (IntToNodeMap::reverse_iterator it2 = neighborhood_sizes2.rbegin(), end2 = neighborhood_sizes2.rend(); it2 != end2; ++it2) {
      Node v = it2->second;
      if (reduced_nodes[u][_mol2.get_id(v)])
        continue;
      if (_mol1.get_color(u) == _mol2.get_color(v)) {
        if (!has_canon[v]) {
          generate_subgraph_canonization(_mol2, v, subgraph_canons);
          has_canon[v] = true;
        }
        Canonization& canon2 = subgraph_canons[v];
        if (are_isomorphic(canon1, canon2)) {
          // generate product node uv
          const Node& uv = add_node(u, v);
          // apply neighborhood subset rule to neighbors of u and v
          BitSet &neighbors_u = neighborhoods1[u];
          const NodeVector& order1 = canon1.get_node_order();
          const NodeVector& order2 = canon2.get_node_order();

          ShortToNodeVectorPairMap current_reductions;
          bfs_subgraph(_mol1, uv, u, neighbors_u, neighborhoods1, order1, order2, current_reductions);
          for (ShortToNodeVectorPairMap::const_iterator it3 = current_reductions.begin(), end3 = current_reductions.end(); it3 != end3; ++it3) {
            for (int i = 0; i < it3->second.first.size(); ++i) {
              for (int j = 0; j < it3->second.second.size(); ++j) {
                reduced_nodes[it3->second.first[i]][_mol2.get_id(it3->second.second[j])] = true;
              }
            }
          }
        }
      }
    }
  }

}

void Product::generate_edges() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
    Node u1 = _g_to_mol1[u1v1];
    Node v1 = _g_to_mol2[u1v1];
    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2) {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _g_to_mol1[u2v2];
      Node v2 = _g_to_mol2[u2v2];

      assert(_mol1.get_color(u1) == _mol2.get_color(v1));

      if (u1 != u2 && v1 != v2) {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 == v1v2) {
          _connectivity[_g.addEdge(u1v1, u2v2)] = u1u2;
        }
      }
    }
  }

}

Node Product::add_node(const Node &u, const Node &v) {
  Node uv = _g.addNode();
  _g_to_mol1[uv] = u;
  _g_to_mol2[uv] = v;
  return uv;
}

//void Product::get_node_mapping(const NodeVectorVector &cliques, NodePairVectorVector &mapping) {
//  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
//    NodePairVector pairs;
//    for (NodeVector::const_iterator it2 = it->begin(), end2 = it->end(); it2 != end2; ++it2) {
//      Node u = _g_to_mol1[*it2];
//      Node v = _g_to_mol2[*it2];
//      pairs.push_back(std::make_pair(u,v));
//      for (NodePairVector::const_iterator it3 = _reductions[*it2].begin(), end3 = _reductions[*it2].end();
//           it3 != end3; ++it3) {
//        pairs.push_back(*it3);
//      }
//    }
//    mapping.push_back(pairs);
//  }
//}
//
//void Product::prune_cliques(NodeVectorVector &cliques) {
//
//  Orbits orbits1(_mol1);
//  Orbits orbits2(_mol2);
//
//  // is comparing the orbit sets enough? We might need to check for induced subgraph isomorphism (of fragments including shells)!
//  IntVectorPairVector orbits;
//  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
//    IntVectorPair pair = std::make_pair<IntVector, IntVector>(IntVector(), IntVector());
//    for (NodeVector::const_iterator it2 = it->begin(), end2 = it->end(); it2 != end2; ++it2) {
//      Node u = _g_to_mol1[*it2];
//      Node v = _g_to_mol2[*it2];
//      pair.first.push_back(orbits1.get_orbit_id(u));
//      pair.second.push_back(orbits2.get_orbit_id(v));
//    }
//    orbits.push_back(pair);
//  }
//
//  for (size_t i = 0; i < cliques.size(); ++i) {
//    for (size_t j = cliques.size() - 1; j >= i; --j) {
//      bool subset = false;
//      bool superset = false;
//      for (size_t k = 0, l = 0; k < cliques[i].size(), l < cliques[j].size();) {
//        if (orbits[i].first[k] == orbits[j].first[l] && orbits[i].second[k] == orbits[j].second[l]) {
//          ++k;
//          ++l;
//        }
//        if (k == cliques[i].size()) {
//          subset = true;
//        } else if (l == cliques[j].size()) {
//          superset = true;
//        }
//      }
//
//      if (subset) {
//        cliques.erase(cliques.begin()+i);
//        --i;
//        break;
//      } else if (superset) {
//        cliques.erase(cliques.begin()+j);
//        break;
//      }
//    }
//  }
//
//}


