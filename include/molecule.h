#include <utility>

#include <utility>

#include <utility>

//
// Created by M. Engler on 10/21/16.
//

#ifndef MOGLI_MOLECULE_H
#define MOGLI_MOLECULE_H

#include <map>
#include <istream>
#include <ostream>
#include <assert.h>
#include <lemon/adaptors.h>
#include <lemon/lgf_reader.h>
#include <lemon/connectivity.h>
#include <queue>
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
        _id_property(std::move(id_property)),
        _color_property(std::move(color_property)),
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

    LGFIOConfig& add_bool_node_prop(const std::string& property) {
      _bool_node_props.insert(property);
      return *this;
    }

    LGFIOConfig& add_int_node_prop(const std::string& property) {
      _int_node_props.insert(property);
      return *this;
    }

    LGFIOConfig& add_double_node_prop(const std::string& property) {
      _double_node_props.insert(property);
      return *this;
    }

    LGFIOConfig& add_string_node_prop(const std::string& property) {
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

  class Molecule {

  public:

  protected:

    typedef typename Graph::template NodeMap<Any> NodeToAnyMap;

    typedef typename UniquePtrMap<std::string, NodeToAnyMap>::type StringToAnyTypeMapMap;

    typedef typename Graph::template NodeMap<unsigned short> NodeToUShortMap;
    typedef typename Graph::template NodeMap<double> NodeToDoubleMap;
    typedef typename Graph::template NodeMap<std::string> NodeToStringMap;

    typedef typename UniquePtrMap<int, Node>::type IntToNodeMap;

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

    explicit Molecule(PeriodicTable &periodic_table) : _g(),
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

    void get_properties(StringVector& properties) const {
      for (const auto & it : _properties) {
        properties.push_back(it.first);
      }

    }

    // set properties
    const Any get_property(Node node, const std::string& property) const {
      return _properties.at(property)->operator[](node);
    }

    void set_property(Node node, const std::string& property, Any value) {
      if (_properties.count(property) == 0) {
        _properties[property] = std::make_unique<NodeToAnyMap>(_g);
      }
      _properties[property]->operator[](node) = std::move(value);
    }

    // add atoms & edges

    const Node add_atom(std::string element) {
      return add_atom(_perdiodic_table.get_number(std::move(element)));
    }

    const Node add_atom(int id, std::string element) {
      return add_atom(id, _perdiodic_table.get_number(std::move(element)));
    }

    const Node add_atom(unsigned short color) {
      return add_atom(_max_uid+1, color);
    }

    const Node add_atom(int id, unsigned short color);

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
      return *_id_to_node.at(id);
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

    void write_gml_stream(const std::string& label, std::ostream &out, bool raw = false);

    std::string write_gml(const std::string& label);

    void write_lgf_stream(std::ostream &out);

    void write_lgf_stream(std::ostream &out, const LGFIOConfig& config);

    std::string write_lgf();

    std::string write_lgf(const LGFIOConfig &config);

    void read_lgf_stream(std::istream &in);

    void read_lgf_stream(std::istream &in, const LGFIOConfig &config);

    void read_lgf(const std::string &in);

    void read_lgf(const std::string &in, const LGFIOConfig &config);

    // connected?

    const bool is_connected() {
      if (_is_connected < 0) {
        _is_connected = is_connected0();
      }
      return _is_connected == 1;
    }

    void get_connected_components(SharedPtrVector<Molecule>::type &components);

    int split(int max_size, int shell, SharedPtrVector<Molecule>::type &components);

    const bool is_isomorphic(Molecule &other) const;

    virtual const std::string print_dot() const;

    const std::string print_dot(const StringVector &properties) const;

    virtual const void print_dot(std::ostream& out) const;

    const void print_dot(std::ostream &out, const StringVector &properties) const;

  protected:

    void balanced_cut(const std::shared_ptr<NodeToBoolMap>& subgraph_map,
                      int max_size, int shell,
                      SharedPtrVector<NodeToBoolMap>::type &partitions);

    int get_bctree(const std::shared_ptr<SubGraph>& subgraph,
                   Graph& bctree,
                   NodeToIntMap& bctree_node_weights, NodeToIntMap& bctree_2_bicon_comp,
                   SharedPtrVector<std::set<Node>>::type & biconnected_nodes,
                   EdgeToIntSetMap& bctree_shell_ids_u, EdgeToIntSetMap& bctree_shell_ids_v, int shell);

    void get_subgraph(Graph& bctree,
                      NodeToIntMap& bctree_components,
                      NodeToIntMap& bctree_2_bicon_comp,
                      SharedPtrVector<std::set<Node>>::type &biconnected_nodes,
                      IntSet& shell_ids,
                      int component,
                      const std::shared_ptr<NodeToBoolMap>& subgraph_map);

    void find_shell_nodes(const std::shared_ptr<SubGraph>& subgraph,
                          NodeToIntMap& bctree_2_bicon_comp,
                          SharedPtrVector<std::set<Node>>::type &biconnected_nodes,
                          Graph& bctree,
                          int shell,
                          EdgeToIntSetMap& bctree_shell_ids_u,
                          EdgeToIntSetMap& bctree_shell_ids_v);

    void bfs_find_shell(const std::shared_ptr<SubGraph>& subgraph,
                        const Node &v, int shell,
                        SubGraph::NodeMap<IntSet>& shell_ids);

    Edge get_cut_edge(Graph& bctree, NodeToIntMap &node_weights,
                      EdgeToIntSetMap& shell_ids_u, EdgeToIntSetMap& shell_ids_v, int N);

    virtual short is_connected0();

    void dfs(const Node& current, NodeToBoolMap& visited);

  };

}

#endif //MOGLI_MOLECULE_H
