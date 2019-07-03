////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
//    Copyright (C) 2016-2019  Martin S. Engler                                                                       //
//                                                                                                                    //
//    This program is free software: you can redistribute it and/or modify                                            //
//    it under the terms of the GNU Lesser General Public License as published                                        //
//    by the Free Software Foundation, either version 3 of the License, or                                            //
//    (at your option) any later version.                                                                             //
//                                                                                                                    //
//    This program is distributed in the hope that it will be useful,                                                 //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                  //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                    //
//    GNU General Public License for more details.                                                                    //
//                                                                                                                    //
//    You should have received a copy of the GNU Lesser General Public License                                        //
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../include/molecule.h"
#include "../include/canonization.h"

using namespace mogli;

// public methods

const Node Molecule::add_atom(int id, unsigned short color)  {
  Node n = _g.addNode();
  _node_to_id[n] = id;
  _id_to_node[id] = std::make_unique<Node>(n);
  if (id > _max_uid) {
    _max_uid = id;
  }
  _atom_count++;
  _colors[n] = color;
  _is_connected = -1;
  return n;
}

std::string Molecule::write_gml(const std::string& label) {
  std::stringstream buffer;
  write_gml_stream(label, buffer, false);
  return buffer.str();
}

std::string Molecule::write_lgf(const LGFIOConfig &config) {
  std::stringstream buffer;
  write_lgf_stream(buffer, config);
  return buffer.str();
}

std::string Molecule::write_lgf() {
  std::stringstream buffer;
  write_lgf_stream(buffer);
  return buffer.str();
}

void Molecule::write_lgf_stream(std::ostream &out) {
  assert(out.good());
  write_lgf_stream(out, LGFIOConfig::get_default());
}

void Molecule::write_lgf_stream(std::ostream &out, const LGFIOConfig &config) {
  assert(out.good());
  // nodes header
  out << "@nodes\n";
  out << config.get_id_property() << "\t" << config.get_color_property() << "\t";
  for (const std::string& prop : config.get_bool_node_props()) {
    out << prop << "\t";
  }
  for (const std::string& prop : config.get_int_node_props()) {
    out << prop << "\t";
  }
  for (const std::string& prop : config.get_double_node_props()) {
    out << prop << "\t";
  }
  for (const std::string& prop : config.get_string_node_props()) {
    out << prop << "\t";
  }
  out << "\n";

  // nodes
  IntVector nodes;
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    int id = _node_to_id[v];
    nodes.push_back(id);
  }
  std::sort(nodes.begin(), nodes.end());

  for (int id : nodes) {
    auto v = *_id_to_node[id];
    out << id << "\t" << _colors[v] << "\t";
    for (const std::string& prop : config.get_bool_node_props()) {
      out << std::get<bool>(get_property(v, prop)) << "\t";
    }
    for (const std::string& prop : config.get_int_node_props()) {
      out << std::get<int>(get_property(v, prop)) << "\t";
    }
    for (const std::string& prop : config.get_double_node_props()) {
      out << std::get<double>(get_property(v, prop)) << "\t";
    }
    for (const std::string& prop : config.get_string_node_props()) {
      out << std::get<std::string>(get_property(v, prop)) << "\t";
    }
    out << "\n";
  }

  // edges
  out << "@edges\n\t\tlabel\t\n";
  std::vector<std::pair<int, int> > edges;
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    int u = _node_to_id[_g.u(e)];
    int v = _node_to_id[_g.v(e)];
    edges.emplace_back(u, v);
  }
  std::sort(edges.begin(), edges.end(), sort_tuple());
  int k = 0;
  for (std::pair<int, int> pair : edges) {
    out << pair.first << "\t" << pair.second << "\t" << k++ << "\t\n";
  }
}

void Molecule::read_lgf(const std::string &in, const LGFIOConfig &config) {
  std::stringstream buffer;
  buffer.str(in);
  read_lgf_stream(buffer, config);
}

void Molecule::read_lgf(const std::string &in) {
  std::stringstream buffer;
  buffer.str(in);
  read_lgf_stream(buffer);
}

void Molecule::read_lgf_stream(std::istream &in) {
  assert(in.good());
  read_lgf_stream(in, LGFIOConfig::get_default());
}

void Molecule::read_lgf_stream(std::istream &in, const LGFIOConfig &config) {
  assert(in.good());

  lemon::GraphReader<Graph> reader(_g, in);
  reader.nodeMap(config.get_color_property(), _colors);
  reader.nodeMap(config.get_id_property(), _node_to_id);

  UniquePtrMap<std::string, NodeToBoolMap>::type bool_props;
  UniquePtrMap<std::string, NodeToIntMap>::type int_props;
  UniquePtrMap<std::string, NodeToDoubleMap>::type double_props;
  UniquePtrMap<std::string, NodeToStringMap>::type string_props;

  for (const std::string& prop : config.get_bool_node_props()) {
    bool_props[prop] = std::make_unique<NodeToBoolMap>(_g);
    reader.nodeMap(prop, *bool_props[prop]);
  }
  for (const std::string& prop : config.get_int_node_props()) {
    int_props[prop] = std::make_unique<NodeToIntMap>(_g);
    reader.nodeMap(prop, *int_props[prop]);
  }
  for (const std::string& prop : config.get_double_node_props()) {
    double_props[prop] = std::make_unique<NodeToDoubleMap>(_g);
    reader.nodeMap(prop, *double_props[prop]);
  }
  for (const std::string& prop : config.get_string_node_props()) {
    string_props[prop] = std::make_unique<NodeToStringMap>(_g);
    reader.nodeMap(prop, *string_props[prop]);
  }
  reader.run();

  _atom_count = 0;
  _max_uid = 0;
  for (NodeIt n(_g); n != lemon::INVALID; ++n) {
    int id = _node_to_id[n];
    _id_to_node[id] = std::make_unique<Node>(n);
    if (id > _max_uid) {
      _max_uid = id;
    }
    _atom_count++;
  }

  for (auto & it : bool_props) {
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      set_property(n, it.first, it.second->operator[](n));
    }
  }
  for (auto & it : int_props) {
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      set_property(n, it.first, it.second->operator[](n));
    }
  }
  for (auto & it : double_props) {
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      set_property(n, it.first, it.second->operator[](n));
    }
  }
  for (auto & it : string_props) {
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      set_property(n, it.first, it.second->operator[](n));
    }
  }

}

void Molecule::write_gml_stream(const std::string& label, std::ostream &out, bool raw) {
  assert(out.good());
  if (!raw) {
    out << "graph [\n\tdirected 0\n\tlabel "<< label << "\n";
  } else {
    out << R"(graph [\n\tdirected 0\n\tlabel )" << label << R"(\n)";
  }


  // nodes
  IntVector nodes;
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    int id = _node_to_id[v];
    nodes.push_back(id);
  }
  std::sort(nodes.begin(), nodes.end());
  if (!raw) {
    for (int id : nodes) {
      out << "\tnode [\n\t\tid " << id << "\n\t\tlabel \"" << _colors[*_id_to_node[id]] << "\"\n\t]\n";
    }
  } else {
    for (int id : nodes) {
      out << R"(\tnode [\n\t\tid )" << id << R"(\n\t\tlabel \")" << _colors[*_id_to_node[id]] << R"(\"\n\t]\n)";
    }
  }

  // edges
  std::vector<std::pair<int, int> > edges;
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    int u = _node_to_id[_g.u(e)];
    int v = _node_to_id[_g.v(e)];
    edges.emplace_back(u, v);
  }
  std::sort(edges.begin(), edges.end(), sort_tuple());
  if (!raw) {
    for (std::pair<int, int> pair : edges) {
      out << "\tedge [\n\t\tsource " << pair.first << "\n\t\ttarget " << pair.second
          << "\n\t\tlabel \"-\"\n\t]\n";
    }
  } else {
    for (std::pair<int, int> pair : edges) {
      out << R"(\tedge [\n\t\tsource )" << pair.first << R"(\n\t\ttarget )" << pair.second
          << R"(\n\t\tlabel \"-\"\n\t]\n)";
    }
  }
  if (!raw) {
    out << "]\n\n";
  } else {
    out << R"(]\n\n)";
  }
}

const bool Molecule::is_isomorphic(Molecule &other) const {
  if (_atom_count != other.get_atom_count())
        return false;
      Canonization c1(*this);
      Canonization c2(other);
      return c1.is_isomorphic(c2);
}

void Molecule::get_connected_components(mogli::SharedPtrVector<mogli::Molecule>::type &components)  {
  NodeToIntMap compMap(_g);
  int count = lemon::connectedComponents(_g, compMap);
  // create a new molecule for every connected component
  for (int i = 0; i < count; ++i) {
    components.push_back(std::make_shared<Molecule>());
  }
  // copy nodes and properties
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    Node cv = components[compMap[v]]->add_atom(_node_to_id[v], _colors[v]);
    for (auto & it : _properties) {
      components[compMap[v]]->set_property(cv, it.first, get_property(v, it.first));
    }
  }
  // copy edges
  lemon::ArcLookUp<Graph> arcLookUp(_g);
  for (int i = 0; i < count; ++i) {
    for (NodeIt u(components[i]->get_graph()); u != lemon::INVALID; ++u) {
      for (NodeIt v = u; v != lemon::INVALID; ++v) {
        if (u == v)
          continue;
        if (arcLookUp(*_id_to_node[components[i]->get_id(u)], *_id_to_node[components[i]->get_id(v)]) != lemon::INVALID) {
          components[i]->add_edge(u,v);
        }
      }
    }
  }

}

int Molecule::split(int max_size, int shell, mogli::SharedPtrVector<mogli::Molecule>::type &components) {

  auto subgraph_map = std::make_shared<NodeToBoolMap>(_g, true);

  SharedPtrVector<NodeToBoolMap>::type partitions;

  balanced_cut(subgraph_map, max_size, shell, partitions);

  for (auto & partition : partitions) {

    auto submol = std::make_shared<Molecule>();
    lemon::ArcLookUp<Graph> arcLookUp(_g);
    // copy core nodes and properties
    for (NodeIt v(_g); v != lemon::INVALID; ++v) {
      if (partition->operator[](v)) {
        Node cv = submol->add_atom(_node_to_id[v], _colors[v]);
        for (auto & it : _properties) {
          submol->set_property(cv, it.first, get_property(v, it.first));
        }
      }
    }

    // copy edges
    for (NodeIt u(submol->get_graph()); u != lemon::INVALID; ++u) {
      for (NodeIt v = u; v != lemon::INVALID; ++v) {
        if (u == v)
          continue;
        Node src_u = *_id_to_node.at(submol->get_id(u));
        Node src_v = *_id_to_node.at(submol->get_id(v));

        if (arcLookUp(src_u, src_v) != lemon::INVALID) {
          submol->add_edge(u,v);
        }
      }
    }

    components.push_back(submol);

  }

  return static_cast<int>(components.size());

}

const std::string Molecule::print_dot() const  {
  std::stringstream buffer;
  print_dot(buffer);
  return buffer.str();
}

const std::string Molecule::print_dot(const mogli::StringVector &properties) const {
  std::stringstream buffer;
  print_dot(buffer, properties);
  return buffer.str();
}

const void Molecule::print_dot(std::ostream &out) const {
  // header
  out << "graph G {" << std::endl
      << "\toverlap=scale" << std::endl;

  // nodes
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    out << "\t" << _g.id(v);
    out << "[style=\"filled\",fillcolor=" << _perdiodic_table.get_color(_colors[v]);
    out << ",label=\"" << _node_to_id[v] << "\"]";
    out << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e)) << std::endl;
  }

  out << "}" << std::endl;
}

const void Molecule::print_dot(std::ostream &out, const mogli::StringVector &properties) const {
  // header
  out << "graph G {" << std::endl
      << "\toverlap=scale" << std::endl;

  // nodes
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    out << "\t" << _g.id(v);
    out << "[style=\"filled\",fillcolor=" << _perdiodic_table.get_color(_colors[v]);
    out << ",label=\"";
    bool first = true;
    for (const auto & prop : properties) {
      if (!first) {
        out << ",";
      } else {
        first = false;
      }
      Any value = get_property(v, prop);
      if (std::holds_alternative<bool>(value)) {
        out << std::get<bool>(value);
      } else if (std::holds_alternative<int>(value)) {
        out << std::get<int>(value);
      } else if (std::holds_alternative<double>(value)) {
        out << std::get<double>(value);
      } else if (std::holds_alternative<std::string>(value)) {
        out << std::get<std::string>(value);
      }
    }
    out << "\"]";
    out << std::endl;
  }

  // edges
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e)) << std::endl;
  }

  out << "}" << std::endl;
}

// protected methods

short Molecule::is_connected0() {
  if (_atom_count < 2)
    return 1;

  NodeToBoolMap visited(_g, false);
  NodeIt v = get_node_iter();
  dfs(v, visited);

  for (; v != lemon::INVALID; ++v) {
    if (!visited[v])
      return 0;
  }
  return 1;
}

void Molecule::dfs(const mogli::Node &current, mogli::NodeToBoolMap &visited) {
  visited[current] = true;
  for (IncEdgeIt e = get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
    Node w = get_opposite_node(current, e);
    if (!visited[w]) {
      dfs(w, visited);
    }
  }
}

void Molecule::balanced_cut(const std::shared_ptr<NodeToBoolMap>& subgraph_map,
                            int max_size, int shell,
                            SharedPtrVector<NodeToBoolMap>::type &partitions) {

  EdgeToBoolMap subgraph_edge_map(_g, true);
  auto subgraph = std::make_shared<SubGraph>(_g, *subgraph_map, subgraph_edge_map);

  int size = lemon::countNodes(*subgraph);
  if (size <= max_size) {
    partitions.push_back(subgraph_map);
    return;
  }

  Graph bctree;
  NodeToIntMap bctree_2_bicon_comp(bctree), bctree_node_weights(bctree);
  EdgeToIntSetMap bctree_shell_ids_u(bctree), bctree_shell_ids_v(bctree);
  SharedPtrVector<std::set<Node>>::type biconnected_nodes;

  int bc_count = get_bctree(subgraph, bctree,
                            bctree_node_weights, bctree_2_bicon_comp, biconnected_nodes,
                            bctree_shell_ids_u, bctree_shell_ids_v, shell);
  if (bc_count == 0) {
    partitions.push_back(subgraph_map);
    return;
  }

  int num_blocks = 0;
  for (NodeIt v(bctree); v != lemon::INVALID; ++v) {
    if (bctree_2_bicon_comp[v] >= 0) {
      ++num_blocks;
    }
  }

  if (num_blocks < 2) {
    partitions.push_back(subgraph_map);
    return;
  }

  assert(lemon::connected(bctree));

  Edge cut_edge  = get_cut_edge(bctree, bctree_node_weights, bctree_shell_ids_u, bctree_shell_ids_v, size);
  if (cut_edge == lemon::INVALID) {
    partitions.push_back(subgraph_map);
    return;
  }

  NodeToBoolMap stnm(bctree, true);
  EdgeToBoolMap stem(bctree, true);
  stem[cut_edge] = false;
  SubGraph subtree(bctree, stnm, stem);
  SubGraph::NodeMap<int> st_components(subtree);
  int count = lemon::connectedComponents(subtree, st_components);

  assert(count == 2);

  auto subgraph_map_0 = std::make_shared<NodeToBoolMap>(_g, false);
  auto subgraph_map_1 = std::make_shared<NodeToBoolMap>(_g, false);

  if (st_components[bctree.u(cut_edge)] == 0) {
    get_subgraph(bctree, st_components,
                 bctree_2_bicon_comp, biconnected_nodes, bctree_shell_ids_u[cut_edge],
                 0, subgraph_map_0);
    get_subgraph(bctree, st_components,
                 bctree_2_bicon_comp, biconnected_nodes, bctree_shell_ids_v[cut_edge],
                 1, subgraph_map_1);
  } else {
    get_subgraph(bctree, st_components,
                 bctree_2_bicon_comp, biconnected_nodes, bctree_shell_ids_v[cut_edge],
                 0, subgraph_map_0);
    get_subgraph(bctree, st_components,
                 bctree_2_bicon_comp, biconnected_nodes, bctree_shell_ids_u[cut_edge],
                 1, subgraph_map_1);
  }

  balanced_cut(subgraph_map_0, max_size, shell, partitions);
  balanced_cut(subgraph_map_1, max_size, shell, partitions);
}

int Molecule::get_bctree(const std::shared_ptr<SubGraph>& subgraph,
                         Graph& bctree,
                         NodeToIntMap& bctree_node_weights, NodeToIntMap& bctree_2_bicon_comp,
                         SharedPtrVector<std::set<Node>>::type & biconnected_nodes,
                         EdgeToIntSetMap& bctree_shell_ids_u, EdgeToIntSetMap& bctree_shell_ids_v, int shell) {
  // get the biconnected components of the graph
  SubGraph::EdgeMap<int> biConnected(*subgraph);
  int bc_count = lemon::biNodeConnectedComponents(*subgraph, biConnected);
  if (bc_count == 0) {
    return 0;
  }

  // get the node sets of the biconnected components
  for (int i = 0; i < bc_count; ++i) {
    biconnected_nodes.push_back(std::make_shared<std::set<Node>> ());
  }
  for (SubGraph::EdgeIt e(*subgraph); e != lemon::INVALID; ++e) {
    int bc_num = biConnected[e];
    (*biconnected_nodes.at(bc_num)).insert(subgraph->u(e));
    (*biconnected_nodes.at(bc_num)).insert(subgraph->v(e));
  }

  // build bc tree
  IntToNodeMap comp_2_bc;
  IntToNodeMap cut_2_bc;
  IntSet cut_ids;

  for (int i = 0; i < biconnected_nodes.size(); ++i) {
    Node u = bctree.addNode();
    bctree_2_bicon_comp[u] = i;
    comp_2_bc[i] = std::make_unique<Node>(u);
    bctree_node_weights[u] = static_cast<int>((*biconnected_nodes[i]).size());
  }

  lemon::ArcLookUp<Graph> arcsBC(bctree);
  for (int i = 0; i < biconnected_nodes.size()-1; ++i) {
    for (int j = i+1; j < biconnected_nodes.size(); ++j) {
      std::vector<Node> intersection;
      std::set<Node> &c1 = *biconnected_nodes.at(i);
      std::set<Node> &c2 = *biconnected_nodes.at(j);

      std::set_intersection(c1.begin(), c1.end(),
                            c2.begin(), c2.end(),
                            std::back_inserter(intersection));

      if (!intersection.empty()) {
        int cut_id = _node_to_id[intersection.front()];
        if (cut_ids.find(cut_id) == cut_ids.end()) {
          Node cut = bctree.addNode();
          bctree_2_bicon_comp[cut] = -_node_to_id[intersection.front()]-1;
          bctree_node_weights[cut] = 1;
          cut_2_bc[cut_id] = std::make_unique<Node>(cut);
          cut_ids.insert(cut_id);
          if (arcsBC(*comp_2_bc[i], cut) == lemon::INVALID) {
            bctree.addEdge(*comp_2_bc[i], cut);
            arcsBC.refresh(*comp_2_bc[i]);
          }
          if (arcsBC(*comp_2_bc[j], cut) == lemon::INVALID) {
            bctree.addEdge(*comp_2_bc[j], cut);
            arcsBC.refresh(*comp_2_bc[j]);
          }
        } else {
          Node cut = *cut_2_bc[cut_id];
          if (arcsBC(*comp_2_bc[i], cut) == lemon::INVALID) {
            bctree.addEdge(*comp_2_bc[i], cut);
            arcsBC.refresh(*comp_2_bc[i]);
          }
          if (arcsBC(*comp_2_bc[j], cut) == lemon::INVALID) {
            bctree.addEdge(*comp_2_bc[j], cut);
            arcsBC.refresh(*comp_2_bc[j]);
          }
        }
      }

    }
  }

  find_shell_nodes(subgraph, bctree_2_bicon_comp, biconnected_nodes, bctree, shell, bctree_shell_ids_u, bctree_shell_ids_v);
  return bc_count;

}

void Molecule::get_subgraph(Graph& bctree,
                            NodeToIntMap& bctree_components,
                            NodeToIntMap& bctree_2_bicon_comp,
                            SharedPtrVector<std::set<Node>>::type &biconnected_nodes,
                            IntSet& shell_ids,
                            int component,
                            const std::shared_ptr<NodeToBoolMap>& subgraph_map) {

  // make union of nodes represented by block nodes in the bc-tree
  for (NodeIt bc_v(bctree); bc_v != lemon::INVALID; ++bc_v) {
    int bc_comp = bctree_2_bicon_comp[bc_v];

    if (bctree_components[bc_v] == component) {
      if (bc_comp >= 0) {
        for (auto & v : *biconnected_nodes[bc_comp]) {
          if (!subgraph_map->operator[](v)) {
            subgraph_map->operator[](v) = true;
          }
        }
      }
    }
  }

  for (int id : shell_ids) {
    auto v = *_id_to_node[id];
    if (!subgraph_map->operator[](v)) {
      subgraph_map->operator[](v) = true;
    }
  }

}

void Molecule::find_shell_nodes(const std::shared_ptr<SubGraph>& subgraph,
                                NodeToIntMap& bctree_2_bicon_comp,
                                SharedPtrVector<std::set<Node>>::type &biconnected_nodes,
                                Graph& bctree,
                                int shell,
                                EdgeToIntSetMap& bctree_shell_ids_u,
                                EdgeToIntSetMap& bctree_shell_ids_v) {

  SubGraph::NodeMap<IntSet> shell_ids(*subgraph);
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    bfs_find_shell(subgraph, v, shell, shell_ids);
  }

  NodeToBoolMap bctree_node_map(bctree, true);
  EdgeToBoolMap bctree_edge_map(bctree, true);
  SubGraph subtree(bctree, bctree_node_map, bctree_edge_map);

  for (EdgeIt e(bctree); e != lemon::INVALID; ++e) {
    bctree_edge_map[e] = false;

    SubGraph::NodeMap<int> bctree_components(subtree);
    int count = lemon::connectedComponents(subtree, bctree_components);
    assert(count == 2);

    IntSet core_nodes_0, core_nodes_1, shell_nodes_0, shell_nodes_1;

    for (SubGraph::NodeIt bc_v(subtree); bc_v != lemon::INVALID; ++bc_v) {
      int bc_comp = bctree_2_bicon_comp[bc_v];

      if (bc_comp >= 0) {
        if (bctree_components[bc_v] == 0) {
          for (auto & v : *biconnected_nodes[bc_comp]) {
            core_nodes_0.insert(_node_to_id[v]);
            for (int sv : shell_ids[v]) {
              shell_nodes_0.insert(sv);
            }
          }
        } else {
          for (auto & v : *biconnected_nodes[bc_comp]) {
            core_nodes_1.insert(_node_to_id[v]);
            for (int sv : shell_ids[v]) {
              shell_nodes_1.insert(sv);
            }
          }
        }
      }
    }

    std::vector<int> diff_0, diff_1;

    std::set_difference(shell_nodes_0.begin(), shell_nodes_0.end(),
                        core_nodes_0.begin(), core_nodes_0.end(),
                        std::back_inserter(diff_0));
    std::set_difference(shell_nodes_1.begin(), shell_nodes_1.end(),
                        core_nodes_1.begin(), core_nodes_1.end(),
                        std::back_inserter(diff_1));

    Node bc_u = bctree.u(e);
    if (bctree_components[bc_u] == 0) {
      for (int id : diff_0) {
        bctree_shell_ids_u[e].insert(id);
      }
      for (int id : diff_1) {
        bctree_shell_ids_v[e].insert(id);
      }

    } else {
      for (int id : diff_1) {
        bctree_shell_ids_u[e].insert(id);
      }
      for (int id : diff_0) {
        bctree_shell_ids_v[e].insert(id);
      }
    }

    bctree_edge_map[e] = true;
  }
}

void Molecule::bfs_find_shell(const std::shared_ptr<SubGraph>& subgraph,
                              const Node &v, int shell,
                              SubGraph::NodeMap<IntSet>& shell_ids) {
  SubGraph::NodeMap<int> depth(*subgraph, 0);
  SubGraph::NodeMap<bool> visited(*subgraph, false);
  std::deque<Node> queue;

  queue.push_back(v);
  visited[v] = true;
  while (!queue.empty()) {
    Node &current = queue.front();
    // if current node is not the core node
    if (current != v) {
      // have we seen the node before?
      int id = _node_to_id[current];
      if (shell_ids[v].count(id) == 0) {
        shell_ids[v].insert(id);
      }
    }

    // breadth-first-search
    if (depth[current] < shell) {
      for (SubGraph::IncEdgeIt e(*subgraph, current); e != lemon::INVALID; ++e) {
        Node w = subgraph->oppositeNode(current, e);
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

Edge Molecule::get_cut_edge(Graph& bctree, NodeToIntMap &node_weights,
                            EdgeToIntSetMap& shell_ids_u, EdgeToIntSetMap& shell_ids_v, int N) {
  NodeToIntMap weight_sum(bctree);
  NodeToIntMap unmarked(bctree);
  EdgeToBoolMap marked(bctree);

  int min_dist = std::numeric_limits<int>::max();
  Edge cut_edge = lemon::INVALID;

  UniquePtrMap<int, std::set<Node>>::type unmarked2node;

  for (NodeIt v(bctree); v != lemon::INVALID; ++v) {
    int deg = 0;
    for (IncEdgeIt e(bctree, v); e != lemon::INVALID; ++e) {
      ++deg;
    }

    if (unmarked2node.find(deg) == unmarked2node.end()) {
      unmarked2node[deg] = std::make_unique<std::set<Node>> ();
    }
    unmarked2node[deg]->insert(v);
    unmarked[v] = deg;
    weight_sum[v] = 0;
  }

  for (EdgeIt e(bctree); e != lemon::INVALID; ++e) {
    marked[e] = false;
  }

  while (!unmarked2node[1]->empty()) {

    std::vector<Node> copy;
    std::copy(unmarked2node[1]->begin(),
              unmarked2node[1]->end(),
              std::back_inserter(copy));

    for (Node u : copy) {
      unmarked2node[1]->erase(u);
      for (IncEdgeIt e(bctree, u); e != lemon::INVALID; ++e) {
        if (!marked[e]) {
          int weight_u, weight_v;
          if (bctree.u(e) == u) {
            weight_u = node_weights[u] + weight_sum[u];
            weight_v = N-weight_u+1;
          } else {
            weight_v = node_weights[u] + weight_sum[u];
            weight_u = N-weight_v+1;
          }

          int wu = weight_u + static_cast<int>(shell_ids_u[e].size());
          int wv = weight_v + static_cast<int>(shell_ids_v[e].size());
          int dist = std::abs(wu - wv);
          if (wu < N && wv < N && dist < min_dist) {
            min_dist = dist;
            cut_edge = e;
          }

          marked[e] = true;

          Node v = bctree.oppositeNode(u, e);
          int num_unmarked = unmarked[v];
          unmarked2node[num_unmarked]->erase(v);

          unmarked[v] = --num_unmarked;
          if (unmarked2node.find(num_unmarked) == unmarked2node.end()) {
            unmarked2node[num_unmarked] = std::make_unique<std::set<Node>> ();
          }
          unmarked2node[num_unmarked]->insert(v);

          if (bctree.u(e) == u) {
            weight_sum[v] = weight_sum[v] + weight_u - 1;
          } else {
            weight_sum[v] = weight_sum[v] + weight_v - 1;
          }

          break;
        }
      }

    }
  }

  return cut_edge;

}
