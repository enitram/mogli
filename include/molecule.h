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
      // FIXME boost::python looses error msg
      try {
        return (&_properties.at(property))->operator[](node);
      } catch (boost::bad_ptr_container_operation) {
        std::string msg = "Could not find property " + property;
        BOOST_PTR_CONTAINER_THROW_EXCEPTION(true, boost::bad_ptr_container_operation, msg.c_str());
      }
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

    virtual void write_gml_stream(std::string label, std::ostream &out);

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
          for (NodeIt v(components[i]->get_graph()); v != lemon::INVALID; ++v) {
            if (u == v)
              continue;
            if (arcLookUp(_id_to_node[components[i]->get_id(u)], _id_to_node[components[i]->get_id(v)]) != lemon::INVALID) {
              components[i]->add_edge(u,v);
            }
          }
        }
      }

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

  inline std::string Molecule::write_gml(std::string label) {
    std::stringstream buffer;
    write_gml_stream(label, buffer);
    return buffer.str();
  }

  inline void Molecule::write_gml_stream(std::string label, std::ostream &out) {
    assert(out.good());
    out << "graph [\n\tdirected 0\n\tlabel "<< label <<"\n";
    for (NodeIt v(_g); v != lemon::INVALID; ++v) {
      out << "\tnode [\n\t\tid " << _node_to_id[v] << "\n\t\tlabel \"" << _colors[v] << "\"\n\t]\n";
    }
    for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
      out << "\tedge [\n\t\tsource " << _node_to_id[_g.u(e)] << "\n\t\ttarget " << _node_to_id[_g.v(e)]
          << "\n\t\tlabel \"-\"\n\t]\n";
    }
    out << "]\n\n";
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
    int k = 0;
    std::vector<std::pair<int, int> > edges;
    for (EdgeIt e(_g); e != lemon::INVALID; ++e, ++k) {
      int u = _node_to_id[_g.u(e)];
      int v = _node_to_id[_g.v(e)];
      edges.push_back(std::make_pair(u, v));
    }
    std::sort(edges.begin(), edges.end(), sort_tuple());
    for (std::pair<int, int> pair : edges) {
      out << pair.first << "\t" << pair.second << "\t" << k << "\t\n";
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
