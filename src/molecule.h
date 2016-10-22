//
// Created by martin on 10/21/16.
//

#ifndef MOGLI_MOLECULE_H
#define MOGLI_MOLECULE_H

#include <istream>
#include <lemon/list_graph.h>
#include <map>
#include <assert.h>
#include <lemon/lgf_reader.h>

namespace mogli {

  typedef lemon::ListGraph Graph;
  typedef Graph::Node Node;
  typedef Graph::NodeIt NodeIt;
  typedef Graph::Edge Edge;
  typedef Graph::EdgeIt EdgeIt;
  typedef Graph::IncEdgeIt IncEdgeIt;

  class Molecule {

  private:

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

  public:
    Molecule() : _g(),
                 _atom_count(0),
                 _label_to_node(),
                 _node_to_label(_g),
                 _node_to_color(_g),
                 _node_to_partial_charge(_g) {}

    ~Molecule() {}

    const Graph &getGraph() const {
      return _g;
    }

    const unsigned int getAtomCount() const {
      return _atom_count;
    }

    const Node getNodeByID(unsigned short id) const {
      return _g.nodeFromId(id);
    }

    const Node getNodeByLabel(const std::string &label) const {
      return _label_to_node.at(label);
    }

    const std::string &getLabel(const Node &node) const {
      return _node_to_label[node];
    }

    const unsigned short getID(const Node &node) const {
      return static_cast<unsigned short>(_g.id(node));
    }

    const unsigned short getColor(const Node &node) const {
      return _node_to_color[node];
    }

    const double getPartialCharge(Node &node) const {
      return _node_to_partial_charge[node];
    }

    virtual void readLGFStream(std::istream &in);

    virtual void readLGF(const std::string &in);

  };

  inline void Molecule::readLGF(const std::string &in) {
    std::stringstream buffer;
    buffer.str(in);
    readLGFStream(buffer);
  }

  inline void Molecule::readLGFStream(std::istream &in) {
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

}

#endif //MOGLI_MOLECULE_H
