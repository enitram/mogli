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
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/any.hpp>
#include "periodictable.h"
#include "types.h"

namespace mogli {

  typedef typename Graph::template NodeMap<int> NodeToIntMap;

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


  class Molecule {

  public:

  protected:

    typedef typename Graph::template NodeMap<boost::any> NodeToAnyMap;

    typedef typename boost::ptr_map<std::string, NodeToAnyMap> StringToAnyTypeMapMap;

    typedef typename Graph::template NodeMap<unsigned short> NodeToUShortMap;
    typedef typename Graph::template NodeMap<double> NodeToDoubleMap;
    typedef typename Graph::template NodeMap<std::string> NodeToStringMap;

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

    // lgf reading

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

//  TODO LGF writer

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
