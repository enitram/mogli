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
#include "iacm.h"

namespace mogli {

  typedef lemon::ListGraph Graph;
  typedef Graph::Node Node;
  typedef Graph::NodeIt NodeIt;
  typedef Graph::Edge Edge;
  typedef Graph::EdgeIt EdgeIt;
  typedef Graph::IncEdgeIt IncEdgeIt;
  typedef std::vector<Node> NodeVector;
  typedef std::vector<std::string> StringVector;
  typedef typename Graph::template NodeMap<bool> NodeToBoolMap;

  class Molecule {

  public:

  private:

    typedef typename Graph::template NodeMap<unsigned short> NodeToUShortMap;
    typedef typename Graph::template NodeMap<int> NodeToIntMap;
    typedef typename Graph::template NodeMap<double> NodeToDoubleMap;
    typedef typename Graph::template NodeMap<std::string> NodeToStringMap;

    typedef typename std::map<std::string, Node> StringToNodeMap;

    typedef typename boost::ptr_map<std::string, NodeToBoolMap> BoolPropertiesMap;
    typedef typename boost::ptr_map<std::string, NodeToIntMap> IntPropertiesMap;
    typedef typename boost::ptr_map<std::string, NodeToDoubleMap> DoublePropertiesMap;
    typedef typename boost::ptr_map<std::string, NodeToStringMap> StringPropertiesMap;

    typedef typename std::map<std::string, StringToNodeMap> InverseStringPropertiesMap;

    Graph _g;

    unsigned int _atom_count;

    NodeToUShortMap _colors;

    BoolPropertiesMap _bool_prop;
    IntPropertiesMap _int_prop;
    DoublePropertiesMap _double_prop;
    StringPropertiesMap _string_prop;

    InverseStringPropertiesMap _inv_string_prop;

    IACM &_iacm;

    short _is_connected;

  public:
    Molecule() : _g(),
                 _colors(_g),
                 _atom_count(0),
                 _is_connected(-1),
                 _iacm(IACM::get_default()){}

    // register properties

    Molecule& add_bool_property(std::string property) {
      _bool_prop.insert(property, new NodeToBoolMap(_g));
      return *this;
    }

    Molecule& add_int_property(std::string property) {
      _int_prop.insert(property, new NodeToIntMap(_g));
      return *this;
    }

    Molecule& add_double_property(std::string property) {
      _double_prop.insert(property, new NodeToDoubleMap(_g));
      return *this;
    }

    Molecule& add_string_property(std::string property) {
      _string_prop.insert(property, new NodeToStringMap(_g));
      _inv_string_prop[property] = StringToNodeMap();
      return *this;
    }

    // get property keys

    const void get_bool_properties(StringVector& keys) const {
      keys.clear();
      keys.reserve(_bool_prop.size());
      for (BoolPropertiesMap::const_iterator it = _bool_prop.begin(), end = _bool_prop.end(); it != end; ++it) {
        keys.push_back(it->first);
      }
    }

    const void get_int_properties(StringVector& keys) const {
      keys.clear();
      keys.reserve(_int_prop.size());
      for (IntPropertiesMap::const_iterator it = _int_prop.begin(), end = _int_prop.end(); it != end; ++it) {
        keys.push_back(it->first);
      }
    }

    const void get_double_properties(StringVector& keys) const {
      keys.clear();
      keys.reserve(_double_prop.size());
      for (DoublePropertiesMap::const_iterator it = _double_prop.begin(), end = _double_prop.end(); it != end; ++it) {
        keys.push_back(it->first);
      }
    }

    const void get_string_properties(StringVector& keys) const {
      keys.clear();
      keys.reserve(_string_prop.size());
      for (StringPropertiesMap::const_iterator it = _string_prop.begin(), end = _string_prop.end(); it != end; ++it) {
        keys.push_back(it->first);
      }
    }

    // set property for node

    const void set_property(Node node, std::string property, bool value) {
      (&_bool_prop.at(property))->set(node, value);
    }

    const void set_property(Node node, std::string property, int value) {
      (&_int_prop.at(property))->set(node, value);
    }

    const void set_property(Node node, std::string property, double value) {
      (&_double_prop.at(property))->set(node, value);
    }

    const void set_property(Node node, std::string property, char* value) {
      set_property(node, property, std::string(value));
    }

    const void set_property(Node node, std::string property, std::string value) {
      (&_string_prop.at(property))->set(node, value);
      if (_inv_string_prop.find(property) != _inv_string_prop.end())
        _inv_string_prop[property][value] = node;
    }

    // get property for node

    const bool get_bool_property(Node node, std::string property) const {
      return (&_bool_prop.at(property))->operator[](node);
    }

    const int get_int_property(Node node, std::string property) const {
      return (&_int_prop.at(property))->operator[](node);
    }

    const double get_double_property(Node node, std::string property) const {
      return (&_double_prop.at(property))->operator[](node);
    }

    const std::string get_string_property(Node node, std::string property) const {
      return (&_string_prop.at(property))->operator[](node);
    }

    // get node for string property

    const Node get_node_by_string_property(std::string property, std::string value) {
      return _inv_string_prop[property][value];
    }

    // add atoms & edges

    const Node add_atom(std::string element) {
      return add_atom(_iacm.get_number(element));
    }

    const Node add_atom(unsigned short color) {
      Node n = _g.addNode();
      _atom_count++;
      _colors[n] = color;
      _is_connected = -1;
      return n;
    }

    const Edge add_edge(Node &u, Node &v) {
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
      return _g.nodeFromId(id);
    }

    const int get_id(const Node &node) const {
      return _g.id(node);
    }

    const unsigned short get_color(const Node &node) const {
      return _colors[node];
    }

    const std::string get_iacm_element(const Node &node) const {
      return _iacm.get_iacm_element(get_color(node));
    }

    const std::string get_chem_element(const Node &node) const {
      return _iacm.get_chem_element(get_color(node));
    }

    const std::string get_chem_color(const Node &node) const {
      return _iacm.get_chem_color(get_color(node));
    }

    // lgf reading

    virtual void read_lgf_stream(std::istream &in);

    virtual void read_lgf(const std::string &in);

    // connected?

    const bool is_connected() {
      if (_is_connected < 0) {
        _is_connected = is_connected0();
      }
      return _is_connected == 1;
    }

    const void print_dot(std::ostream& out, const StringVector& properties = {}) const {
      // header
      out << "graph G {" << std::endl
          << "\toverlap=scale" << std::endl
          << "\tlayout=neato" << std::endl;

      StringVector bool_props, int_props, double_props, string_props;
      if (properties.size() > 0) {
        StringVector keys;
        get_bool_properties(keys);
        check_properties(properties, keys, bool_props);
        get_int_properties(keys);
        check_properties(properties, keys, int_props);
        get_double_properties(keys);
        check_properties(properties, keys, double_props);
        get_string_properties(keys);
        check_properties(properties, keys, string_props);
      }

      // nodes
      for (NodeIt v(_g); v != lemon::INVALID; ++v) {
        out << "\t" << _g.id(v);
        if (properties.size() > 0) {
          out << "[style=filled,fillcolor=" << _iacm.get_chem_color(_colors[v]);
          out << ",label=\"";
          bool first = true;
          for (std::vector<std::string>::const_iterator it = string_props.begin(), end = string_props.end(); it != end; ++it) {
            if (first) {
              out  << get_string_property(v, *it);
              first = false;
            } else {
              out << "\\n" << "," << get_string_property(v, *it);
            }
          }
          for (std::vector<std::string>::const_iterator it = double_props.begin(), end = double_props.end(); it != end; ++it) {
            if (first) {
              out  << get_double_property(v, *it);
              first = false;
            } else {
              out << "\\n" << "," << get_double_property(v, *it);
            }
          }
          for (std::vector<std::string>::const_iterator it = int_props.begin(), end = int_props.end(); it != end; ++it) {
            if (first) {
              out  << get_int_property(v, *it);
              first = false;
            } else {
              out << "\\n" << "," << get_int_property(v, *it);
            }
          }
          for (std::vector<std::string>::const_iterator it = bool_props.begin(), end = bool_props.end(); it != end; ++it) {
            if (first) {
              out  << get_bool_property(v, *it);
              first = false;
            } else {
              out << "\\n" << "," << get_bool_property(v, *it);
            }
          }
          out << "\"]";
        }
        out << std::endl;
      }

      // edges
      for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
        out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e)) << std::endl;
      }

      out << "}" << std::endl;
    }

  private:

    const void inline check_properties(const StringVector& from, const StringVector& keys, StringVector& to) const {
      for (std::vector<std::string>::const_iterator it = from.begin(), end = from.end(); it != end; ++it) {
        if (std::find(keys.begin(), keys.end(), *it) != keys.end()) {
          to.push_back(*it);
        }
      }
    }

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

  inline void Molecule::read_lgf(const std::string &in) {
    std::stringstream buffer;
    buffer.str(in);
    read_lgf_stream(buffer);
  }

  inline void Molecule::read_lgf_stream(std::istream &in) {
    assert(in.good());
    _g.clear();

    lemon::GraphReader<Graph> reader(_g, in);
    reader.nodeMap("atomType", _colors);
    for (BoolPropertiesMap::iterator it = _bool_prop.begin(), end = _bool_prop.end(); it != end; ++it) {
      reader.nodeMap(it->first, *(it->second));
    }
    for (IntPropertiesMap::iterator it = _int_prop.begin(), end = _int_prop.end(); it != end; ++it) {
      reader.nodeMap(it->first, *(it->second));
    }
    for (DoublePropertiesMap::iterator it = _double_prop.begin(), end = _double_prop.end(); it != end; ++it) {
      reader.nodeMap(it->first, *(it->second));
    }
    for (StringPropertiesMap::iterator it = _string_prop.begin(), end = _string_prop.end(); it != end; ++it) {
      reader.nodeMap(it->first, *(it->second));
    }
    reader.run();

    _atom_count = 0;
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      for (StringPropertiesMap::const_iterator it = _string_prop.begin(), end = _string_prop.end(); it != end; ++it) {
        _inv_string_prop[it->first][(&_string_prop.at(it->first))->operator[](n)] = n;
      }
      _atom_count++;
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
