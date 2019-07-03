//
// Created by M. Engler on 26/10/16.
//

#include "../include/product.h"

using namespace mogli;

// constructors

Product::Product(const mogli::Product &parent, int component)
    : _mol1(parent._mol1)
    , _mol2(parent._mol2)
    , _shell(parent._shell)
    , _gen_type(parent._gen_type)
    , _g()
    , _node_sizes(_g)
    , _reductions(_g)
    , _g_to_mol1(_g)
    , _g_to_mol2(_g)
    , _g_to_mol1_canons(_g)
    , _g_to_mol2_canons(_g)
    , _connectivity(_g)
    , _comp_map(_g)
    , _comp_count(0)
    , _comp_sizes()
    , _is_complete(false) {

  NodeToNodeMap nodes(_g);
  int V = 0;
  for (NodeIt orig(parent._g); orig != lemon::INVALID; ++orig) {
    if (parent._comp_map[orig] == component) {
      Node copy = _g.addNode();
      ++V;
      nodes[copy] = orig;
      _node_sizes[copy] = parent._node_sizes[orig];
      _reductions[copy] = parent._reductions[orig];
      _g_to_mol1[copy] = parent._g_to_mol1[orig];
      _g_to_mol2[copy] = parent._g_to_mol2[orig];
      _g_to_mol1_canons[copy] = parent._g_to_mol1_canons[orig];
      _g_to_mol2_canons[copy] = parent._g_to_mol2_canons[orig];
    }
  }

  lemon::ArcLookUp<Graph> arcLookUp(parent._g);
  int E = 0;
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    for (NodeIt u = v; u != lemon::INVALID; ++u) {
      if (v == u)
        continue;

      Edge e = arcLookUp(nodes[u], nodes[v]);
      if (e != lemon::INVALID) {
        _connectivity[_g.addEdge(u,v)] = parent._connectivity[e];
        ++E;
      }
    }
  }
  _is_complete = E == ((V*(V-1))/2);

}

Product::Product(const mogli::Molecule &mol1, const mogli::Molecule &mol2, int shell,
                 mogli::Product::GenerationType gen, unsigned int min_core_size)
    : _mol1(mol1)
    , _mol2(mol2)
    , _shell(shell)
    , _gen_type(gen)
    , _g()
    , _node_sizes(_g)
    , _reductions(_g)
    , _g_to_mol1(_g)
    , _g_to_mol2(_g)
    , _g_to_mol1_canons(_g)
    , _g_to_mol2_canons(_g)
    , _connectivity(_g)
    , _comp_map(_g)
    , _comp_count(1)
    , _is_complete(false) {
  int V, E;
  if ((gen == DEG_1 || gen == UNCON_DEG_1) && shell > 0) {
    V = generate_nodes_deg1();
  } else if ((gen == SUB || gen == UNCON_SUB) && shell > 0) {
    V = generate_nodes_sub();
  } else {
    V = generate_nodes();
  }
  if (gen == UNCON || gen == UNCON_DEG_1 || gen == UNCON_SUB) {
    E = generate_edges_connected(min_core_size);
  } else {
    E = generate_edges();
  }
  _is_complete = E == ((V*(V-1))/2);

}

// public methods

const std::string Product::print_dot() const {
  std::stringstream buffer;
  print_dot(buffer);
  return buffer.str();
}

const std::string Product::print_dot(const mogli::StringVector &properties) const {
  std::stringstream buffer;
  print_dot(buffer, properties);
  return buffer.str();
}

const void Product::print_dot(std::ostream &out) const {
  // header
  out << "graph G {" << std::endl
      << "\toverlap=scale" << std::endl;

  // nodes
  for (NodeIt uv(_g); uv != lemon::INVALID; ++uv) {
    Node u = _g_to_mol1[uv];
    Node v = _g_to_mol2[uv];

    out << "\t" << _g.id(uv);
    out << "[style=\"filled\",fillcolor=" << _mol1.get_color_name(u);
    out << ",label=\"" << _mol1.get_id(u) << "x" << _mol2.get_id(v) << "\"]";
    out << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    out << "\t" << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
    if (_connectivity[e]) {
      out <<  std::endl;
    } else {
      out << "[style=\"dashed\"]" << std::endl;
    }
  }

  out << "}" << std::endl;
}

const void Product::print_dot(std::ostream &out, const mogli::StringVector &properties) const {
  // header
  out << "graph G {" << std::endl
      << "\toverlap=scale" << std::endl;

  // nodes
  for (NodeIt uv(_g); uv != lemon::INVALID; ++uv) {
    Node u = _g_to_mol1[uv];
    Node v = _g_to_mol2[uv];

    out << "\t" << _g.id(uv);
    out << "[style=\"filled\",fillcolor=" << _mol1.get_color_name(u);
    out << ",label=\"";
    bool first = true;
    for (const auto prop : properties) {
      if (!first) {
        out << ",";
      } else {
        first = false;
      }
      Any value1 = _mol1.get_property(u, prop);
      if (std::holds_alternative<bool>(value1)) {
        out << std::get<bool>(value1);
      } else if (std::holds_alternative<int>(value1)) {
        out << std::get<int>(value1);
      } else if (std::holds_alternative<double>(value1)) {
        out << std::get<double>(value1);
      } else if (std::holds_alternative<std::string>(value1)) {
        out << std::get<std::string>(value1);
      }
      out << "x";
      Any value2 = _mol2.get_property(v, prop);
      if (std::holds_alternative<bool>(value2)) {
        out << std::get<bool>(value2);
      } else if (std::holds_alternative<int>(value2)) {
        out << std::get<int>(value2);
      } else if (std::holds_alternative<double>(value2)) {
        out << std::get<double>(value2);
      } else if (std::holds_alternative<std::string>(value2)) {
        out << std::get<std::string>(value2);
      }
    }
    out << "\"]";
    out << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    out << "\t" << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
    if (_connectivity[e]) {
      out <<  std::endl;
    } else {
      out << "[style=\"dashed\"]" << std::endl;
    }
  }

  out << "}" << std::endl;
}

// private methods

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
  while (!queue.empty()) {
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
  while (!queue.empty()) {
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
  while (!queue.empty()) {
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
        ++_node_sizes[product_node];
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

int Product::generate_nodes() {
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
  int V = 0;
  for (NodeIt u = _mol1.get_node_iter(); u != lemon::INVALID; ++u) {
    for (NodeIt v = _mol2.get_node_iter(); v != lemon::INVALID; ++v) {
      if (_mol1.get_color(u) == _mol2.get_color(v) &&
          subgraph_canons1[u].is_isomorphic(subgraph_canons2[v])) {
        const Node& uv = add_node(u, v);
        ++V;
        _node_sizes[uv] = 1;
        _g_to_mol1_canons[uv] = subgraph_canons1[u];
        _g_to_mol2_canons[uv] = subgraph_canons2[v];
      }
    }
  }

  return V;
}

int Product::generate_nodes_deg1() {
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
  int V = 0;
  for (auto it = deg_to_node1.rbegin(), end = deg_to_node1.rend(); it != end; ++it) {
    Node u = it->second;
    NodeToBoolMap filter(g1, false);
    bfs(_mol1, u, filter);
    Canonization canon1 = Canonization(_mol1, filter, u);
    // iterate nodes v with decreasing degree
    for (auto it2 = deg_to_node2.rbegin(), end2 = deg_to_node2.rend(); it2 != end2; ++it2) {
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
          ++V;
          _node_sizes[uv] = 1;
          // apply degree-1 rule to neighbors of u and v
          const ShortVector& order1 = canon1.get_node_order();
          const ShortVector& order2 = canon2.get_node_order();
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
          for (const auto & it3 : current_reductions) {
            NodeVector foo = it3.second.first;
            NodeVector bar = it3.second.second;
            for (int i = 0; i < it3.second.first.size(); ++i) {
              for (int j = 0; j < it3.second.second.size(); ++j) {
                Node _u = it3.second.first[i];
                Node _v = it3.second.second[j];
                reduced_nodes[_u].insert(_mol2.get_id(_v));
              }
            }
          }
        }
      }
    }
  }

  return V;

}

int Product::generate_nodes_sub() {
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
  int V = 0;
  for (auto it = neighborhood_sizes1.rbegin(), end = neighborhood_sizes1.rend(); it != end; ++it) {
    Node u = it->second;
    NodeToBoolMap filter(g1, false);
    bfs(_mol1, u, filter);
    Canonization canon1 = Canonization(_mol1, filter, u);
    // iterate nodes v with decreasing neighborhood size
    for (auto it2 = neighborhood_sizes2.rbegin(), end2 = neighborhood_sizes2.rend(); it2 != end2; ++it2) {
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
          ++V;
          _node_sizes[uv] = 1;
          // apply neighborhood subset rule to neighbors of u and v
          IntSet &neighbors_u = neighborhoods1[u];
          const ShortVector& _order1 = canon1.get_node_order();
          const ShortVector& _order2 = canon2.get_node_order();
          _g_to_mol1_canons[uv] = canon1;
          _g_to_mol2_canons[uv] = canon2;

          NodeVector order1, order2;

          for (const auto & it3 : _order1) {
            order1.push_back(_mol1.get_node_by_id(it3));
          }
          for (const auto it3 : _order2) {
            order2.push_back(_mol2.get_node_by_id(it3));
          }

          ShortToNodeVectorPairMap current_reductions;
          bfs_subgraph(_mol1, uv, u, neighbors_u, neighborhoods1, order1, order2, current_reductions);
          for (const auto & it3 : current_reductions) {
            for (int i = 0; i < it3.second.first.size(); ++i) {
              for (int j = 0; j < it3.second.second.size(); ++j) {
                reduced_nodes[it3.second.first[i]].insert(_mol2.get_id(it3.second.second[j]));
              }
            }
          }
        }
      }
    }
  }

  return V;
}

int Product::generate_edges() {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  int E = 0;
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
    Node u1 = _g_to_mol1[u1v1];
    Node v1 = _g_to_mol2[u1v1];
    assert(_mol1.get_color(u1) == _mol2.get_color(v1));

    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2) {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _g_to_mol1[u2v2];
      Node v2 = _g_to_mol2[u2v2];

      assert(_mol1.get_color(u2) == _mol2.get_color(v2));

      if (u1 != u2 && v1 != v2) {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 == v1v2) {
          _connectivity[_g.addEdge(u1v1, u2v2)] = u1u2;
          ++E;
        }
      }
    }
  }

  return E;
}

int Product::generate_edges_connected(unsigned int min_core_size) {
  const Graph& g1 = _mol1.get_graph();
  const Graph& g2 = _mol2.get_graph();

  lemon::ArcLookUp<Graph> arcLookUp1(g1);
  lemon::ArcLookUp<Graph> arcLookUp2(g2);

  // build c-edges
  int E = 0;
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
    Node u1 = _g_to_mol1[u1v1];
    Node v1 = _g_to_mol2[u1v1];
    assert(_mol1.get_color(u1) == _mol2.get_color(v1));

    for (NodeIt u2v2 = u1v1; u2v2 != lemon::INVALID; ++ u2v2) {
      if (u1v1 == u2v2)
        continue;

      Node u2 = _g_to_mol1[u2v2];
      Node v2 = _g_to_mol2[u2v2];
      assert(_mol1.get_color(u2) == _mol2.get_color(v2));

      if (u1 != u2 && v1 != v2) {
        bool u1u2 = arcLookUp1(u1, u2) != lemon::INVALID;
        bool v1v2 = arcLookUp2(v1, v2) != lemon::INVALID;

        if (u1u2 && v1v2) {
          _connectivity[_g.addEdge(u1v1, u2v2)] = true;
          ++E;
        }
      }
    }
  }

  // find connected components
  NodeToIntMap compMap(_g);
  int count = lemon::connectedComponents(_g, compMap);

  for (int c = 0; c < count; ++c) {

    int comp_size = 0;
    // get the nodes of the current component
    NodeVector nodes;
    for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
      if (compMap[u1v1] != c)
        continue;

      comp_size += _node_sizes[u1v1];
      nodes.push_back(u1v1);
    }

    // delete current component if it is too small
    if (comp_size < min_core_size) {
      for (auto node : nodes) {
        _g.erase(node);
      }
      continue;
    }

    // add d-edges for current component
    for (int i = 0; i < nodes.size()-1; ++i) {
      Node u1v1 = nodes[i];
      Node u1 = _g_to_mol1[u1v1];
      Node v1 = _g_to_mol2[u1v1];
      assert(_mol1.get_color(u1) == _mol2.get_color(v1));

      for (int j = i+1; j < nodes.size(); ++j) {
        Node u2v2 = nodes[j];

        Node u2 = _g_to_mol1[u2v2];
        Node v2 = _g_to_mol2[u2v2];
        assert(_mol1.get_color(u2) == _mol2.get_color(v2));

        if (u1 != u2 && v1 != v2) {
          bool u1u2 = arcLookUp1(u1, u2) == lemon::INVALID;
          bool v1v2 = arcLookUp2(v1, v2) == lemon::INVALID;

          if (u1u2 && v1v2) {
            _connectivity[_g.addEdge(u1v1, u2v2)] = false;
            ++E;
          }
        }
      }
    }

  }

  // count connected components again (we might have deleted some components)
  _comp_count = lemon::connectedComponents(_g, _comp_map);
  for (int c = 0; c < _comp_count; ++c) {
    _comp_sizes.push_back(0);
  }
  for (NodeIt u1v1(_g); u1v1 != lemon::INVALID; ++u1v1) {
    _comp_sizes[_comp_map[u1v1]] += _node_sizes[u1v1];
  }

  return E;
}

Node Product::add_node(const Node &u, const Node &v) {
  Node uv = _g.addNode();
  _g_to_mol1[uv] = u;
  _g_to_mol2[uv] = v;
  return uv;
}



