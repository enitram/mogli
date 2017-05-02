//
// Created by M. Engler on 26/10/16.
//

#include "../include/product.h"

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

void Product::bfs_neighbors(const Molecule &mol, const Node &v, IntSet &neighbors, int& size) {
  NodeToBoolMap visited(mol.get_graph(), false);
  NodeToIntMap depth(mol.get_graph(), 0);
  NodeDeque queue;

  queue.push_back(v);
  visited[v] = true;
  while (queue.size() > 0) {
    Node &current = queue.front();
    neighbors.insert(mol.get_id(current));
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

void Product::bfs_subgraph(const Molecule &mol, const Node &product_node, const Node &root_node, const IntSet &root_neighbors,
                           const NodeToIntSetMap &neighborhoods, const NodeVector &order1, const NodeVector &order2,
                           ShortToNodeVectorPairMap &current_reductions) {
  NodeToBoolMap visited(mol.get_graph(), false);
  NodeToIntMap depth(mol.get_graph(), 0);
  NodeDeque queue;

  queue.push_back(root_node);
  visited[root_node] = true;
  while (queue.size() > 0) {
    Node &current = queue.front();

    const IntSet &neighbors_current = neighborhoods[current];
    if (std::includes(root_neighbors.begin(), root_neighbors.end(),
                      neighbors_current.begin(), neighbors_current.end())) {
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

void Product::generate_subgraph(const Molecule &mol, const Node &v, NodeToIntSetMap &neighborhoods, IntToNodeMap &sizes) {
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
          subgraph_canons1[u].is_isomorphic(subgraph_canons2[v])) {
        const Node& uv = add_node(u, v);
        _node_sizes[uv] = 1;
        _g_to_mol1_canons[uv] = subgraph_canons1[u];
        _g_to_mol2_canons[uv] = subgraph_canons2[v];
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

  NodeToIntSetMap reduced_nodes(g1, IntSet());

  // iterate nodes u with decreasing degree
  for (IntToNodeMap::reverse_iterator it = deg_to_node1.rbegin(), end = deg_to_node1.rend(); it != end; ++it) {
    Node u = it->second;
    NodeToBoolMap filter(g1, false);
    bfs(_mol1, u, filter);
    Canonization canon1 = Canonization(_mol1, filter, u);
    // iterate nodes v with decreasing degree
    for (IntToNodeMap::reverse_iterator it2 = deg_to_node2.rbegin(), end2 = deg_to_node2.rend(); it2 != end2; ++it2) {
      Node v = it2->second;
      if (reduced_nodes[u].count(_mol2.get_id(v)) > 0) {
        continue;
      }
      if (_mol1.get_color(u) == _mol2.get_color(v)) {
        if (!has_canon[v]) {
          generate_subgraph_canonization(_mol2, v, subgraph_canons);
          has_canon[v] = true;
        }
        Canonization& canon2 = subgraph_canons[v];
        if (canon1.is_isomorphic(canon2)) {
          // generate product node uv
          const Node& uv = add_node(u, v);
          _node_sizes[uv] = 1;
          // apply degree-1 rule to neighbors of u and v
          const ShortVector& order1 = canon1.get_node_order();
          const ShortVector& order2 = canon2.get_node_order();
          int u_id = _mol1.get_id(u);
          int v_id = _mol2.get_id(v);
          _g_to_mol1_canons[uv] = canon1;
          _g_to_mol2_canons[uv] = canon2;

          ShortToNodeVectorPairMap current_reductions;
          for (IncEdgeIt e(g1, u); e != lemon::INVALID; ++e) {
            Node w = g1.oppositeNode(u, e);
            if (deg1[w] == 1) {
              int index = static_cast<int>(std::find(order1.begin(), order1.end(), _mol1.get_id(w)) - order1.begin());
              Node x = _mol2.get_node_by_id(order2[index]);
              unsigned short color = _mol1.get_color(w);
              if (current_reductions.find(color) == current_reductions.end()) {
                current_reductions[color] = std::make_pair(NodeVector(), NodeVector());
              }
              current_reductions[color].first.push_back(w);
              current_reductions[color].second.push_back(x);
              _reductions[uv].push_back(std::make_pair(w, x));
              ++_node_sizes[uv];
            }
          }
          for (ShortToNodeVectorPairMap::const_iterator it3 = current_reductions.begin(), end3 = current_reductions.end(); it3 != end3; ++it3) {
            NodeVector foo = it3->second.first;
            NodeVector bar = it3->second.second;
            for (int i = 0; i < it3->second.first.size(); ++i) {
              for (int j = 0; j < it3->second.second.size(); ++j) {
                Node u = it3->second.first[i];
                Node v = it3->second.second[j];
                reduced_nodes[it3->second.first[i]].insert(_mol2.get_id(it3->second.second[j]));
              }
            }
          }
        }
      }
    }
  }

}

// FIXME missing: nodes sizes
void Product::generate_nodes_sub() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  IntToNodeMap neighborhood_sizes1;
  IntToNodeMap neighborhood_sizes2;

  NodeToIntSetMap neighborhoods1(g1, IntSet());
  NodeToIntSetMap neighborhoods2(g2, IntSet());

  for (NodeIt v(g1); v != lemon::INVALID; ++v) {
    generate_subgraph(_mol1, v, neighborhoods1, neighborhood_sizes1);
  }

  for (NodeIt v(g2); v != lemon::INVALID; ++v) {
    generate_subgraph(_mol2, v, neighborhoods2, neighborhood_sizes2);
  }

  NodeToCanonizationMap subgraph_canons(g2);
  NodeToBoolMap has_canon(g2, false);

  NodeToIntSetMap reduced_nodes(g1, IntSet());

  // iterate nodes u with decreasing neighborhood size
  for (IntToNodeMap::reverse_iterator it = neighborhood_sizes1.rbegin(), end = neighborhood_sizes1.rend(); it != end; ++it) {
    Node u = it->second;
    NodeToBoolMap filter(g1, false);
    bfs(_mol1, u, filter);
    Canonization canon1 = Canonization(_mol1, filter, u);
    // iterate nodes v with decreasing neighborhood size
    for (IntToNodeMap::reverse_iterator it2 = neighborhood_sizes2.rbegin(), end2 = neighborhood_sizes2.rend(); it2 != end2; ++it2) {
      Node v = it2->second;
      if (reduced_nodes[u].count(_mol2.get_id(v)) > 0) {
        continue;
      }
      if (_mol1.get_color(u) == _mol2.get_color(v)) {
        if (!has_canon[v]) {
          generate_subgraph_canonization(_mol2, v, subgraph_canons);
          has_canon[v] = true;
        }
        Canonization& canon2 = subgraph_canons[v];
        if (canon1.is_isomorphic(canon2)) {
          // generate product node uv
          const Node& uv = add_node(u, v);
          _node_sizes[uv] = 1;
          // apply neighborhood subset rule to neighbors of u and v
          IntSet &neighbors_u = neighborhoods1[u];
          const ShortVector& _order1 = canon1.get_node_order();
          const ShortVector& _order2 = canon2.get_node_order();
          _g_to_mol1_canons[uv] = canon1;
          _g_to_mol2_canons[uv] = canon2;

          NodeVector order1, order2;

          for (ShortVector::const_iterator it3 = _order1.begin(), end3 = _order1.end(); it3 != end3; ++it3) {
            order1.push_back(_mol1.get_node_by_id(*it3));
          }
          for (ShortVector::const_iterator it3 = _order2.begin(), end3 = _order2.end(); it3 != end3; ++it3) {
            order2.push_back(_mol2.get_node_by_id(*it3));
          }

          ShortToNodeVectorPairMap current_reductions;
          bfs_subgraph(_mol1, uv, u, neighbors_u, neighborhoods1, order1, order2, current_reductions);
          for (ShortToNodeVectorPairMap::const_iterator it3 = current_reductions.begin(), end3 = current_reductions.end(); it3 != end3; ++it3) {
            for (int i = 0; i < it3->second.first.size(); ++i) {
              for (int j = 0; j < it3->second.second.size(); ++j) {
                reduced_nodes[it3->second.first[i]].insert(_mol2.get_id(it3->second.second[j]));
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

void Product::generate_edges(NodeToBoolMap& connected_nodes) {
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
          if (u1u2) {
            connected_nodes[u1v1] = true;
            connected_nodes[u2v2] = true;
          }
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

void Product::prune_nodes(NodeToBoolMap& connected_nodes, unsigned int min_core_size, unsigned int max_core_size) {
  for (NodeIt u(_g); u != lemon::INVALID;) {
    Node current = u;
    ++u;
    if (_gen_type == UNCON_DEG_1 && !connected_nodes[current]) {
      if (_reductions[current].size() == 0 ||
          _node_sizes[current] < min_core_size || _node_sizes[current] > max_core_size) {
        // TODO mark erased nodes/edges instead of real erase?
        _g.erase(current);
      }
    }
  }
}


