//
// Created by M. Engler on 10/21/16.
//

#ifndef MOGLI_MOLECULE_H
#define MOGLI_MOLECULE_H

#include <map>
#include <istream>
#include <ostream>
#include <assert.h>
#include <lemon/core.h>
#include <lemon/adaptors.h>
#include <lemon/lgf_reader.h>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/ptr_container/exception.hpp>
#include <boost/any.hpp>
#include <boost/smart_ptr/make_shared.hpp>
#include <queue>
#include <boost/ptr_container/ptr_vector.hpp>
#include "periodictable.h"
#include "types.h"

namespace mogli {

  typedef typename Graph::template NodeMap<int> NodeToIntMap;
  typedef typename Graph::template EdgeMap<bool> EdgeToBoolMap;
  typedef typename Graph::template EdgeMap<IntSet> EdgeToIntSetMap;
  typedef lemon::SubGraph<Graph, NodeToBoolMap, EdgeToBoolMap> SubGraph;

  class LGFIOConfig {

  public:

    LGFIOConfig(std::string id_property, std::string color_property) :
        _id_property(id_property),
        _color_property(color_property),
        _bool_node_props(),
        _int_node_props(),
        _double_node_props(),
        _string_node_props() {}

  private:

    std::string _id_property;
    std::string _color_property;
    StringSet _bool_node_props;
    StringSet _int_node_props;
    StringSet _double_node_props;
    StringSet _string_node_props;

    LGFIOConfig() :
      _id_property(std::string("label")),
      _color_property(std::string("atomType")),
      _bool_node_props(),
      _int_node_props(),
      _double_node_props(),
      _string_node_props() {
      add_double_node_prop("partial_charge");
      add_string_node_prop("label2");
      add_double_node_prop("coordX");
      add_double_node_prop("coordY");
      add_double_node_prop("coordZ");
      add_int_node_prop("initColor");
    }

  public:

    LGFIOConfig& add_bool_node_prop(std::string property) {
      _bool_node_props.insert(property);
      return *this;
    }

    LGFIOConfig& add_int_node_prop(std::string property) {
      _int_node_props.insert(property);
      return *this;
    }

    LGFIOConfig& add_double_node_prop(std::string property) {
      _double_node_props.insert(property);
      return *this;
    }

    LGFIOConfig& add_string_node_prop(std::string property) {
      _string_node_props.insert(property);
      return *this;
    }

    const std::string &get_id_property() const {
      return _id_property;
    }

    const std::string &get_color_property() const {
      return _color_property;
    }

    const StringSet &get_bool_node_props() const {
      return _bool_node_props;
    }

    const StringSet &get_int_node_props() const {
      return _int_node_props;
    }

    const StringSet &get_double_node_props() const {
      return _double_node_props;
    }

    const StringSet &get_string_node_props() const {
      return _string_node_props;
    }

    static LGFIOConfig &get_default() {
      static LGFIOConfig instance;
      return instance;
    }

  };

  struct sort_tuple {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
      if (left.first == right.first) {
        return left.second < right.second;
      }
      return left.first < right.first;
    }
  };

  template <typename T1, typename T2>
  struct more_first {
    bool operator() (const std::pair<T1, T2>& x,
                     const std::pair<T1, T2>& y) const {
      return x.first > y.first;
    }
  };

  class Molecule {

  public:

  protected:

    typedef typename Graph::template NodeMap<boost::any> NodeToAnyMap;

    typedef typename boost::ptr_map<std::string, NodeToAnyMap> StringToAnyTypeMapMap;

    typedef typename Graph::template NodeMap<unsigned short> NodeToUShortMap;
    typedef typename Graph::template NodeMap<double> NodeToDoubleMap;
    typedef typename Graph::template NodeMap<std::string> NodeToStringMap;
    typedef typename Graph::template NodeMap<char *> NodeToCStringMap;

    typedef typename boost::ptr_map<int, Node> IntToNodeMap;

    Graph _g;

    unsigned int _atom_count;

    NodeToUShortMap _colors;
    NodeToIntMap _node_to_id;
    IntToNodeMap _id_to_node;
    int _max_uid;

    StringToAnyTypeMapMap _properties;

    PeriodicTable &_perdiodic_table;

    short _is_connected;

  public:
    Molecule() : _g(),
                 _colors(_g),
                 _node_to_id(_g),
                 _max_uid(0),
                 _atom_count(0),
                 _is_connected(-1),
                 _properties(),
                 _perdiodic_table(PeriodicTable::get_default()) {}

    Molecule(PeriodicTable &periodic_table) : _g(),
                                             _colors(_g),
                                             _node_to_id(_g),
                                             _max_uid(0),
                                             _atom_count(0),
                                             _is_connected(-1),
                                             _properties(),
                                             _perdiodic_table(periodic_table) {}

    PeriodicTable &get_perdiodic_table() const {
      return _perdiodic_table;
    }

    const void get_properties(StringVector& properties) const {
      for (StringToAnyTypeMapMap::const_iterator it = _properties.begin(), end = _properties.end(); it != end; ++it) {
        properties.push_back(it->first);
      }

    }

    // set properties
    const boost::any get_property(Node node, std::string property) const {
      return (&_properties.at(property))->operator[](node);
    }

    void set_property(Node node, std::string property, boost::any value) {
      if (_properties.count(property) == 0) {
        _properties.insert(property, new NodeToAnyMap(_g));
      }
      (&_properties.at(property))->set(node, value);
    }

    // add atoms & edges

    const Node add_atom(std::string element) {
      return add_atom(_perdiodic_table.get_number(element));
    }

    const Node add_atom(int id, std::string element) {
      return add_atom(id, _perdiodic_table.get_number(element));
    }

    const Node add_atom(unsigned short color) {
      return add_atom(_max_uid+1, color);
    }

    const Node add_atom(int id, unsigned short color) {
      Node n = _g.addNode();
      _node_to_id[n] = id;
      _id_to_node[id] = n;
      if (id > _max_uid) {
        _max_uid = id;
      }
      _atom_count++;
      _colors[n] = color;
      _is_connected = -1;
      return n;
    }

    const Edge add_edge(const Node &u, const Node &v) {
      _is_connected = -1;
      return _g.addEdge(u,v);
    }

    // getters & iterators

    const unsigned int get_atom_count() const {
      return _atom_count;
    }

    const Graph &get_graph() const {
      return _g;
    }

    const NodeIt get_node_iter() const {
      return NodeIt(_g);
    }

    const EdgeIt get_edge_iter() const {
      return EdgeIt(_g);
    }

    const IncEdgeIt get_inc_edge_iter(const Node &node) const {
      return IncEdgeIt(_g, node);
    }

    const Node get_opposite_node(const Node &node, const Edge &edge) const {
      return _g.oppositeNode(node, edge);
    }

    const Node get_u(const Edge &edge) const {
      return _g.u(edge);
    }

    const Node get_v(const Edge &edge) const {
      return _g.v(edge);
    }

    const bool has_node_with_id(int id) const {
      return _id_to_node.find(id) != _id_to_node.end();
    }

    const Node get_node_by_id(int id) const {
      return _id_to_node.at(id);
    }

    const int get_id(const Node &node) const {
      return _node_to_id[node];
    }

    const unsigned short get_color(const Node &node) const {
      return _colors[node];
    }

    const std::string get_element(const Node &node) const {
      return _perdiodic_table.get_element(get_color(node));
    }

    const std::string get_color_name(const Node &node) const {
      return _perdiodic_table.get_color(get_color(node));
    }

    // I/O

    virtual void write_gml_stream(std::string label, std::ostream &out, bool raw = false);

    virtual std::string write_gml(std::string label);

    virtual void write_lgf_stream(std::ostream &out);

    virtual void write_lgf_stream(std::ostream &out, const LGFIOConfig& config);

    virtual std::string write_lgf();

    virtual std::string write_lgf(const LGFIOConfig &config);

    virtual void read_lgf_stream(std::istream &in);

    virtual void read_lgf_stream(std::istream &in, const LGFIOConfig &config);

    virtual void read_lgf(const std::string &in);

    virtual void read_lgf(const std::string &in, const LGFIOConfig &config);

    // connected?

    const bool is_connected() {
      if (_is_connected < 0) {
        _is_connected = is_connected0();
      }
      return _is_connected == 1;
    }

    void get_connected_components(std::vector<boost::shared_ptr<Molecule> > &components) {
      NodeToIntMap compMap(_g);
      int count = lemon::connectedComponents(_g, compMap);
      // create a new molecule for every connected component
      for (int i = 0; i < count; ++i) {
        boost::shared_ptr<Molecule> mol = boost::make_shared<Molecule>();
        components.push_back(mol);
      }
      // copy nodes and properties
      for (NodeIt v(_g); v != lemon::INVALID; ++v) {
        Node cv = components[compMap[v]]->add_atom(_node_to_id[v], _colors[v]);
        for (StringToAnyTypeMapMap::const_iterator it = _properties.begin(), end = _properties.end(); it != end; ++it) {
          components[compMap[v]]->set_property(cv, it->first, get_property(v, it->first));
        }
      }
      // copy edges
      lemon::ArcLookUp<Graph> arcLookUp(_g);
      for (int i = 0; i < count; ++i) {
        for (NodeIt u(components[i]->get_graph()); u != lemon::INVALID; ++u) {
          for (NodeIt v = u; v != lemon::INVALID; ++v) {
            if (u == v)
              continue;
            if (arcLookUp(_id_to_node[components[i]->get_id(u)], _id_to_node[components[i]->get_id(v)]) != lemon::INVALID) {
              components[i]->add_edge(u,v);
            }
          }
        }
      }

    }

    int split(int max_size, int shell, std::vector<boost::shared_ptr<Molecule> > &components) {

      boost::shared_ptr<NodeToBoolMap> subgraph_map = boost::make_shared<NodeToBoolMap>(_g, true);

      std::vector<boost::shared_ptr<NodeToBoolMap> > partitions;

      balanced_cut(subgraph_map, max_size, shell, partitions);

      for (int i = 0; i < partitions.size(); ++i) {

        boost::shared_ptr<Molecule> submol = boost::make_shared<Molecule>();
        lemon::ArcLookUp<Graph> arcLookUp(_g);
        // copy core nodes and properties
        for (NodeIt v(_g); v != lemon::INVALID; ++v) {
          if (partitions.at(i)->operator[](v)) {
            Node cv = submol->add_atom(_node_to_id[v], _colors[v]);
            for (StringToAnyTypeMapMap::const_iterator it = _properties.begin(), end = _properties.end(); it != end; ++it) {
              submol->set_property(cv, it->first, get_property(v, it->first));
            }
          }
        }

        // copy edges
        for (NodeIt u(submol->get_graph()); u != lemon::INVALID; ++u) {
          for (NodeIt v = u; v != lemon::INVALID; ++v) {
            if (u == v)
              continue;
            Node src_u = _id_to_node.at(submol->get_id(u));
            Node src_v = _id_to_node.at(submol->get_id(v));

            if (arcLookUp(src_u, src_v) != lemon::INVALID) {
              submol->add_edge(u,v);
            }
          }
        }

        components.push_back(submol);

      }

      return static_cast<int>(components.size());

    }

  protected:

    void balanced_cut(boost::shared_ptr<NodeToBoolMap> subgraph_map,
                      int max_size, int shell,
                      std::vector<boost::shared_ptr<NodeToBoolMap> > &partitions) {

      EdgeToBoolMap subgraph_edge_map(_g, true);
      boost::shared_ptr<SubGraph> subgraph = boost::make_shared<SubGraph>(_g, *subgraph_map, subgraph_edge_map);

      int size = lemon::countNodes(*subgraph);
      if (size <= max_size) {
        partitions.push_back(subgraph_map);
        return;
      }

      Graph bctree;
      NodeToIntMap bctree_2_bicon_comp(bctree), bctree_node_weights(bctree);
      EdgeToIntSetMap bctree_shell_ids_u(bctree), bctree_shell_ids_v(bctree);
      boost::ptr_vector<std::set<Node> > biconnected_nodes;

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

      boost::shared_ptr<NodeToBoolMap> subgraph_map_0 = boost::make_shared<NodeToBoolMap>(_g, false);
      boost::shared_ptr<NodeToBoolMap> subgraph_map_1 = boost::make_shared<NodeToBoolMap>(_g, false);

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

    int get_bctree(boost::shared_ptr<SubGraph> subgraph,
                    Graph& bctree,
                    NodeToIntMap& bctree_node_weights, NodeToIntMap& bctree_2_bicon_comp,
                    boost::ptr_vector<std::set<Node> >& biconnected_nodes,
                    EdgeToIntSetMap& bctree_shell_ids_u, EdgeToIntSetMap& bctree_shell_ids_v, int shell) {
      // get the biconnected components of the graph
      SubGraph::EdgeMap<int> biConnected(*subgraph);
      int bc_count = lemon::biNodeConnectedComponents(*subgraph, biConnected);
      if (bc_count == 0) {
        return 0;
      }

      // get the node sets of the biconnected components
      for (int i = 0; i < bc_count; ++i) {
        biconnected_nodes.push_back(new std::set<Node> ());
      }
      for (SubGraph::EdgeIt e(*subgraph); e != lemon::INVALID; ++e) {
        int bc_num = biConnected[e];
        biconnected_nodes.at(bc_num).insert(subgraph->u(e));
        biconnected_nodes.at(bc_num).insert(subgraph->v(e));
      }

      // build bc tree
      IntToNodeMap comp_2_bc;
      IntToNodeMap cut_2_bc;
      IntSet cut_ids;

      for (int i = 0; i < biconnected_nodes.size(); ++i) {
        Node u = bctree.addNode();
        bctree_2_bicon_comp[u] = i;
        comp_2_bc[i] = u;
        bctree_node_weights[u] = static_cast<int>(biconnected_nodes.at(i).size());
      }

      lemon::ArcLookUp<Graph> arcsBC(bctree);
      for (int i = 0; i < biconnected_nodes.size()-1; ++i) {
        for (int j = i+1; j < biconnected_nodes.size(); ++j) {
          std::vector<Node> intersection;
          std::set<Node> &c1 = biconnected_nodes.at(i);
          std::set<Node> &c2 = biconnected_nodes.at(j);

          std::set_intersection(c1.begin(), c1.end(),
                                c2.begin(), c2.end(),
                                std::back_inserter(intersection));

          if (intersection.size() > 0) {
            int cut_id = _node_to_id[intersection.front()];
            if (cut_ids.find(cut_id) == cut_ids.end()) {
              Node cut = bctree.addNode();
              bctree_2_bicon_comp[cut] = -_node_to_id[intersection.front()]-1;
              bctree_node_weights[cut] = 1;
              cut_2_bc[cut_id] = cut;
              cut_ids.insert(cut_id);
              if (arcsBC(comp_2_bc[i], cut) == lemon::INVALID) {
                bctree.addEdge(comp_2_bc[i], cut);
                arcsBC.refresh(comp_2_bc[i]);
              }
              if (arcsBC(comp_2_bc[j], cut) == lemon::INVALID) {
                bctree.addEdge(comp_2_bc[j], cut);
                arcsBC.refresh(comp_2_bc[j]);
              }
            } else {
              Node cut = cut_2_bc[cut_id];
              if (arcsBC(comp_2_bc[i], cut) == lemon::INVALID) {
                bctree.addEdge(comp_2_bc[i], cut);
                arcsBC.refresh(comp_2_bc[i]);
              }
              if (arcsBC(comp_2_bc[j], cut) == lemon::INVALID) {
                bctree.addEdge(comp_2_bc[j], cut);
                arcsBC.refresh(comp_2_bc[j]);
              }
            }
          }

        }
      }

      find_shell_nodes(subgraph, bctree_2_bicon_comp, biconnected_nodes, bctree, shell, bctree_shell_ids_u, bctree_shell_ids_v);
      return bc_count;

    }

    void get_subgraph(Graph& bctree,
                     NodeToIntMap& bctree_components,
                     NodeToIntMap& bctree_2_bicon_comp,
                     boost::ptr_vector<std::set<Node> >& biconnected_nodes,
                     IntSet& shell_ids,
                     int component,
                     boost::shared_ptr<NodeToBoolMap> subgraph_map) {

      // make union of nodes represented by block nodes in the bc-tree
      for (NodeIt bc_v(bctree); bc_v != lemon::INVALID; ++bc_v) {
        int bc_comp = bctree_2_bicon_comp[bc_v];

        if (bctree_components[bc_v] == component) {
          if (bc_comp >= 0) {
            for (Node v : biconnected_nodes[bc_comp]) {
              if (!subgraph_map->operator[](v)) {
                subgraph_map->operator[](v) = true;
              }
            }
          }
        }
      }

      for (int id : shell_ids) {
        Node v = _id_to_node[id];
        if (!subgraph_map->operator[](v)) {
          subgraph_map->operator[](v) = true;
        }
      }

    }

    void find_shell_nodes(boost::shared_ptr<SubGraph> subgraph,
                          NodeToIntMap& bctree_2_bicon_comp,
                          boost::ptr_vector<std::set<Node> >& biconnected_nodes,
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
              for (Node v : biconnected_nodes[bc_comp]) {
                core_nodes_0.insert(_node_to_id[v]);
                for (int sv : shell_ids[v]) {
                  shell_nodes_0.insert(sv);
                }
              }
            } else {
              for (Node v : biconnected_nodes[bc_comp]) {
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

    void bfs_find_shell(boost::shared_ptr<SubGraph> subgraph,
                        const Node &v, int shell,
                        SubGraph::NodeMap<IntSet>& shell_ids) {
      SubGraph::NodeMap<int> depth(*subgraph, 0);
      SubGraph::NodeMap<bool> visited(*subgraph, false);
      std::deque<Node> queue;

      queue.push_back(v);
      visited[v] = true;
      while (queue.size() > 0) {
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

    Edge get_cut_edge(Graph& bctree, NodeToIntMap &node_weights,
                      EdgeToIntSetMap& shell_ids_u, EdgeToIntSetMap& shell_ids_v, int N) {
      NodeToIntMap weight_sum(bctree);
      NodeToIntMap unmarked(bctree);
      EdgeToBoolMap marked(bctree);

      int min_dist = std::numeric_limits<int>::max();
      Edge cut_edge = lemon::INVALID;

      boost::ptr_map<int, std::set<Node> > unmarked2node;

      for (NodeIt v(bctree); v != lemon::INVALID; ++v) {
        int deg = 0;
        for (IncEdgeIt e(bctree, v); e != lemon::INVALID; ++e) {
          ++deg;
        }

        if (unmarked2node.find(deg) == unmarked2node.end()) {
          unmarked2node[deg] = std::set<Node> ();
        }
        unmarked2node[deg].insert(v);
        unmarked[v] = deg;
        weight_sum[v] = 0;
      }

      for (EdgeIt e(bctree); e != lemon::INVALID; ++e) {
        marked[e] = false;
      }

      while (unmarked2node[1].size() > 0) {

        std::vector<Node> copy;
        std::copy(unmarked2node[1].begin(),
                  unmarked2node[1].end(),
                  std::back_inserter(copy));

        for (Node u : copy) {
          unmarked2node[1].erase(u);
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
              unmarked2node[num_unmarked].erase(v);

              unmarked[v] = --num_unmarked;
              if (unmarked2node.find(num_unmarked) == unmarked2node.end()) {
                unmarked2node[num_unmarked] = std::set<Node>();
              }
              unmarked2node[num_unmarked].insert(v);

              if (bctree.u(e) == u) {
                weight_sum[v] = weight_sum[v] + weight_u - 1;
              } else {
                weight_sum[v] = weight_sum[v] + weight_v - 1;
              };

              break;
            }
          }

        }
      }

      return cut_edge;

    }

  public:

    const bool is_isomorphic(Molecule &other) const;

    virtual const std::string print_dot() const {
      std::stringstream buffer;
      print_dot(buffer);
      return buffer.str();
    }

    const std::string print_dot(const StringVector &properties) const {
      std::stringstream buffer;
      print_dot(buffer, properties);
      return buffer.str();
    }

    virtual const void print_dot(std::ostream& out) const {
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

    const void print_dot(std::ostream &out, const StringVector &properties) const {
      // header
      out << "graph G {" << std::endl
          << "\toverlap=scale" << std::endl;

      // nodes
      for (NodeIt v(_g); v != lemon::INVALID; ++v) {
        out << "\t" << _g.id(v);
        out << "[style=\"filled\",fillcolor=" << _perdiodic_table.get_color(_colors[v]);
        out << ",label=\"";
        bool first = true;
        for (std::string prop : properties) {
          if (!first) {
            out << ",";
          } else {
            first = false;
          }
          boost::any value = get_property(v, prop);
          if (value.type() == typeid(bool)) {
            out << boost::any_cast<bool>(value);
          } else if (value.type() == typeid(int)) {
            out << boost::any_cast<int>(value);
          } else if (value.type() == typeid(double)) {
            out << boost::any_cast<double>(value);
          } else if (value.type() == typeid(std::string)) {
            out << boost::any_cast<std::string>(value);
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

  protected:

    virtual short is_connected0();

    void dfs(const Node& current, NodeToBoolMap& visited) {
      visited[current] = true;
      for (IncEdgeIt e = get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
        Node w = get_opposite_node(current, e);
        if (!visited[w]) {
          dfs(w, visited);
        }
      }
    }

  };

  inline std::string Molecule::write_gml(std::string label) {
    std::stringstream buffer;
    write_gml_stream(label, buffer, false);
    return buffer.str();
  }

  inline void Molecule::write_gml_stream(std::string label, std::ostream &out, bool raw) {
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
        out << "\tnode [\n\t\tid " << id << "\n\t\tlabel \"" << _colors[_id_to_node[id]] << "\"\n\t]\n";
      }
    } else {
      for (int id : nodes) {
        out << R"(\tnode [\n\t\tid )" << id << R"(\n\t\tlabel \")" << _colors[_id_to_node[id]] << R"(\"\n\t]\n)";
      }
    }

    // edges
    std::vector<std::pair<int, int> > edges;
    for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
      int u = _node_to_id[_g.u(e)];
      int v = _node_to_id[_g.v(e)];
      edges.push_back(std::make_pair(u, v));
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

  inline std::string Molecule::write_lgf(const LGFIOConfig &config) {
    std::stringstream buffer;
    write_lgf_stream(buffer, config);
    return buffer.str();
  }

  inline std::string Molecule::write_lgf() {
    std::stringstream buffer;
    write_lgf_stream(buffer);
    return buffer.str();
  }

  inline void Molecule::write_lgf_stream(std::ostream &out) {
    assert(out.good());
    write_lgf_stream(out, LGFIOConfig::get_default());
  }

  inline void Molecule::write_lgf_stream(std::ostream &out, const LGFIOConfig &config) {
    assert(out.good());
    // nodes header
    out << "@nodes\n";
    out << config.get_id_property() << "\t" << config.get_color_property() << "\t";
    for (std::string prop : config.get_bool_node_props()) {
      out << prop << "\t";
    }
    for (std::string prop : config.get_int_node_props()) {
      out << prop << "\t";
    }
    for (std::string prop : config.get_double_node_props()) {
      out << prop << "\t";
    }
    for (std::string prop : config.get_string_node_props()) {
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
      Node v = _id_to_node[id];
      out << id << "\t" << _colors[v] << "\t";
      for (std::string prop : config.get_bool_node_props()) {
        out << boost::any_cast<bool>(get_property(v, prop)) << "\t";
      }
      for (std::string prop : config.get_int_node_props()) {
        out << boost::any_cast<int>(get_property(v, prop)) << "\t";
      }
      for (std::string prop : config.get_double_node_props()) {
        out << boost::any_cast<double>(get_property(v, prop)) << "\t";
      }
      for (std::string prop : config.get_string_node_props()) {
        out << boost::any_cast<std::string>(get_property(v, prop)) << "\t";
      }
      out << "\n";
    }

    // edges
    out << "@edges\n\t\tlabel\t\n";
    std::vector<std::pair<int, int> > edges;
    for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
      int u = _node_to_id[_g.u(e)];
      int v = _node_to_id[_g.v(e)];
      edges.push_back(std::make_pair(u, v));
    }
    std::sort(edges.begin(), edges.end(), sort_tuple());
    int k = 0;
    for (std::pair<int, int> pair : edges) {
      out << pair.first << "\t" << pair.second << "\t" << k++ << "\t\n";
    }
  }

  inline void Molecule::read_lgf(const std::string &in, const LGFIOConfig &config) {
    std::stringstream buffer;
    buffer.str(in);
    read_lgf_stream(buffer, config);
  }

  inline void Molecule::read_lgf(const std::string &in) {
    std::stringstream buffer;
    buffer.str(in);
    read_lgf_stream(buffer);
  }

  inline void Molecule::read_lgf_stream(std::istream &in) {
    assert(in.good());
    read_lgf_stream(in, LGFIOConfig::get_default());
  }

  inline void Molecule::read_lgf_stream(std::istream &in, const LGFIOConfig &config) {
    assert(in.good());

    lemon::GraphReader<Graph> reader(_g, in);
    reader.nodeMap(config.get_color_property(), _colors);
    reader.nodeMap(config.get_id_property(), _node_to_id);

    boost::ptr_map<std::string, NodeToBoolMap> bool_props;
    boost::ptr_map<std::string, NodeToIntMap> int_props;
    boost::ptr_map<std::string, NodeToDoubleMap> double_props;
    boost::ptr_map<std::string, NodeToStringMap> string_props;

    for (std::string prop : config.get_bool_node_props()) {
      bool_props.insert(prop, new NodeToBoolMap(_g));
      reader.nodeMap(prop, *(&bool_props.at(prop)));
    }
    for (std::string prop : config.get_int_node_props()) {
      int_props.insert(prop, new NodeToIntMap(_g));
      reader.nodeMap(prop, *(&int_props.at(prop)));
    }
    for (std::string prop : config.get_double_node_props()) {
      double_props.insert(prop, new NodeToDoubleMap(_g));
      reader.nodeMap(prop, *(&double_props.at(prop)));
    }
    for (std::string prop : config.get_string_node_props()) {
      string_props.insert(prop, new NodeToStringMap(_g));
      reader.nodeMap(prop, *(&string_props.at(prop)));
    }
    reader.run();

    _atom_count = 0;
    _max_uid = 0;
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      int id = _node_to_id[n];
      _id_to_node[id] = n;
      if (id > _max_uid) {
        _max_uid = id;
      }
      _atom_count++;
    }

    for (boost::ptr_map<std::string, NodeToBoolMap>::iterator it = bool_props.begin(), end = bool_props.end();
         it != end; ++it) {
      for (NodeIt n(_g); n != lemon::INVALID; ++n) {
        set_property(n, it->first, it->second->operator[](n));
      }
    }
    for (boost::ptr_map<std::string, NodeToIntMap>::iterator it = int_props.begin(), end = int_props.end();
         it != end; ++it) {
      for (NodeIt n(_g); n != lemon::INVALID; ++n) {
        set_property(n, it->first, it->second->operator[](n));
      }
    }
    for (boost::ptr_map<std::string, NodeToDoubleMap>::iterator it = double_props.begin(), end = double_props.end();
         it != end; ++it) {
      for (NodeIt n(_g); n != lemon::INVALID; ++n) {
        set_property(n, it->first, it->second->operator[](n));
      }
    }
    for (boost::ptr_map<std::string, NodeToStringMap>::iterator it = string_props.begin(), end = string_props.end();
         it != end; ++it) {
      for (NodeIt n(_g); n != lemon::INVALID; ++n) {
        set_property(n, it->first, it->second->operator[](n));
      }
    }

  }

  inline short Molecule::is_connected0() {
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

}

#endif //MOGLI_MOLECULE_H
