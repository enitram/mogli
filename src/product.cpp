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

void Product::dfs(const Molecule &mol, const Node &v, int depth, NodeToBoolMap &visited, NodeToBoolMap &filter) {
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

void Product::dfs(const Molecule &mol, const Node &v, int depth, NodeToBoolMap &visited, BitSet &neighbors, int &size) {
  visited[v] = true;
  ++size;
  neighbors[mol.get_id(v)] = true;
  if (depth < _shell) {
    for (IncEdgeIt e = mol.get_inc_edge_iter(v); e != lemon::INVALID; ++e) {
      Node w = mol.get_opposite_node(v, e);
      if (!visited[w]) {
        dfs(mol, w, depth+1, visited, neighbors, size);
      }
    }
  }
}

void Product::dfs_sub(const Molecule &mol, const Node &product_node, const Node &v, int depth,
                      const BitSet &root_neighbors,
                      const NodeToBitSetMap &neighborhoods, const NodeVector &order1, const NodeVector &order2,
                      NodeToBitSetMap &reduced_nodes, NodeToBoolMap &visited) {
  visited[v] = true;
  const BitSet &neighbors_v = neighborhoods[v];
  if (neighbors_v.is_subset_of(root_neighbors)) {
    int index = static_cast<int>(std::find(order1.begin(), order1.end(), v) - order1.begin());
    Node w = order2[index];
    reduced_nodes[v][_mol2.get_id(w)] = true;
    _reductions[product_node].push_back(std::make_pair(v, w));
    if (depth < _shell) {
      for (IncEdgeIt e = mol.get_inc_edge_iter(v); e != lemon::INVALID; ++e) {
        Node x = mol.get_opposite_node(v, e);
        if (!visited[x]) {
          dfs_sub(mol, product_node, x, depth + 1, root_neighbors, neighborhoods, order1, order2, reduced_nodes,
                  visited);
        }
      }
    }
  }
}

void Product::generate_subgraph_canonization(const Molecule &mol, const Node &v, NodeToCanonizationMap &map) {
  NodeToBoolMap filter(mol.get_graph(), false);
  NodeToBoolMap visited(mol.get_graph(), false);
  dfs(mol, v, 0, visited, filter);
  map[v] = Canonization(mol, filter, v);
}

void Product::generate_subgraph(const Molecule &mol, const Node &v, NodeToBitSetMap &neighborhoods, IntToNodeMap &sizes) {
  NodeToBoolMap visited(mol.get_graph(), false);
  int size = 0;
  dfs(mol, v, 0, visited, neighborhoods[v], size);
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
    NodeToBoolMap visited(g1, false);
    dfs(_mol1, u, 0, visited, filter);
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

          for (IncEdgeIt e(g1, u); e != lemon::INVALID; ++e) {
            Node w = g1.oppositeNode(u, e);
            if (deg1[w] == 1) {
              int index = static_cast<int>(std::find(order1.begin(), order1.end(), w) - order1.begin());
              Node x = order2[index];
              reduced_nodes[w][_mol2.get_id(x)] = true;
              _reductions[uv].push_back(std::make_pair(w, x));
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
    NodeToBoolMap visited(g1, false);
    dfs(_mol1, u, 0, visited, filter);
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
          for (IncEdgeIt e = _mol1.get_inc_edge_iter(u); e != lemon::INVALID; ++e) {
            Node w = _mol1.get_opposite_node(u, e);
            NodeToBoolMap visited2(_mol1.get_graph(), false);
            visited2[u] = true;
            dfs_sub(_mol1, uv, w, 1, neighbors_u, neighborhoods1, order1, order2, reduced_nodes, visited2);
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
