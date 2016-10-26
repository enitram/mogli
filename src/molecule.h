//
// Created by martin on 10/21/16.
//

#ifndef MOGLI_MOLECULE_H
#define MOGLI_MOLECULE_H

#include <map>
#include <istream>
#include <assert.h>
#include <lemon/core.h>
#include <lemon/lgf_reader.h>
#include <lemon/list_graph.h>
#include "iacm.h"

namespace mogli {

  typedef lemon::ListGraph Graph;
  typedef Graph::Node Node;
  typedef Graph::NodeIt NodeIt;
  typedef Graph::Edge Edge;
  typedef Graph::EdgeIt EdgeIt;
  typedef Graph::IncEdgeIt IncEdgeIt;

  class Molecule {

  public:

    static constexpr double UNDEFINED = std::numeric_limits<double>::infinity();

  private:

    typedef typename Graph::template NodeMap<bool> NodeToBoolMap;
    typedef typename std::map<std::string, Node> StringToNodeMap;
    typedef typename Graph::template NodeMap<std::string> NodeToStringMap;
    typedef typename Graph::template NodeMap<double> NodeToDoubleMap;
    typedef typename Graph::template NodeMap<unsigned short> NodeToShortMap;

    Graph _g;

    unsigned int _atom_count;

    StringToNodeMap _label_to_node;
    NodeToStringMap _node_to_label;
    NodeToShortMap _node_to_color;
    NodeToDoubleMap _node_to_partial_charge;
    IACM &_iacm;

    short _is_connected;

  public:
    Molecule() : _g(),
                 _atom_count(0),
                 _label_to_node(),
                 _node_to_label(_g),
                 _node_to_color(_g),
                 _node_to_partial_charge(_g),
                 _is_connected(-1),
                 _iacm(IACM::get_default()){}

    Molecule(const Molecule &molecule) : _g(),
                                         _atom_count(molecule._atom_count),
                                         _label_to_node(),
                                         _node_to_label(_g),
                                         _node_to_color(_g),
                                         _node_to_partial_charge(_g),
                                         _is_connected(molecule._is_connected),
                                         _iacm(molecule._iacm) {
      lemon::GraphCopy<Graph, Graph> copy(molecule._g, _g);
      Graph::NodeMap<Graph::Node> nodeRef(molecule._g);
      copy.nodeRef(nodeRef);
      copy.nodeMap(molecule._node_to_label, _node_to_label);
      copy.nodeMap(molecule._node_to_color, _node_to_color);
      copy.nodeMap(molecule._node_to_partial_charge, _node_to_partial_charge);
      copy.run();
      for (StringToNodeMap::const_iterator it = molecule._label_to_node.begin(), end = molecule._label_to_node.end();
           it != end; ++it) {
        _label_to_node[it->first] = nodeRef[it->second];
      }
    }

    const Node add_atom(std::string element, std::string label = "", double partial_charge = UNDEFINED) {
      return add_atom(_iacm.get_number(element), label, partial_charge);
    }

    const Node add_atom(unsigned short color, std::string label = "", double partial_charge = UNDEFINED) {
      Node n = _g.addNode();
      _atom_count++;
      _node_to_label[n] = label;
      _node_to_color[n] = color;
      _node_to_partial_charge[n] = partial_charge;
      _label_to_node[_node_to_label[n]] = n;
      _is_connected = -1;
      return n;
    }

    const Edge add_edge(Node &u, Node &v) {
      _is_connected = -1;
      return _g.addEdge(u,v);
    }

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

    const Node get_node_by_id(const unsigned short id) const {
      return _g.nodeFromId(id);
    }

    const Node get_node_by_label(const std::string &label) const {
      return _label_to_node.at(label);
    }

    const std::string get_label(const Node &node) const {
      return _node_to_label[node];
    }

    const unsigned short get_id(const Node &node) const {
      return static_cast<unsigned short>(_g.id(node));
    }

    const unsigned short get_color(const Node &node) const {
      return _node_to_color[node];
    }

    const std::string get_iacm_element(const Node &node) const {
      return _iacm.get_iacm_element(get_color(node));
    }

    const std::string get_chem_element(const Node &node) const {
      return _iacm.get_chem_element(get_color(node));
    }

    const double get_partial_charge(const Node &node) const {
      return _node_to_partial_charge[node];
    }

    virtual void read_lgf_stream(std::istream &in);

    virtual void read_lgf(const std::string &in);

    const bool is_connected() {
      if (_is_connected < 0) {
        _is_connected = is_connected0();
      }
      return _is_connected == 1;
    }

  private:

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

    lemon::graphReader(_g, in)
        .nodeMap("partial_charge", _node_to_partial_charge)
        .nodeMap("label", _node_to_label)
        .nodeMap("atomType", _node_to_color)
        .run();

    _atom_count = 0;
    for (NodeIt n(_g); n != lemon::INVALID; ++n) {
      _label_to_node[_node_to_label[n]] = n;
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
