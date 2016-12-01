//
// Created by M. Engler on 10/20/16.
//

#ifndef MOGLI_CANONIZATION_H
#define MOGLI_CANONIZATION_H

#include "molecule.h"
#include <nauty.h>
#include <malloc.h>
#include <boost/dynamic_bitset.hpp>

namespace mogli {

  typedef std::vector<unsigned short> ShortVector;
  typedef std::vector<unsigned long> LongVector;
  typedef boost::dynamic_bitset<> BitSet;

  class Canonization {

  public:

    Canonization() {}

    Canonization(const Molecule& mol) : _colors(), _canonization(), _node_order() {
      init(mol);
    }

    Canonization(const Molecule& mol, const NodeToBoolMap& filter, const Node& root) :
        _colors(), _canonization(), _node_order() {
      init(mol, filter, root);
    }

    Canonization(const ShortVector &_colors, const LongVector &_canonization, const NodeVector &_node_order) :
        _colors(_colors), _canonization(_canonization), _node_order(_node_order) {}

    const ShortVector &get_colors() const {
      return _colors;
    }

    const LongVector &get_canonization() const {
      return _canonization;
    }

    const NodeVector &get_node_order() const {
      return _node_order;
    }

  private:

    typedef std::set<unsigned short> ShortSet;
    typedef typename Graph::template NodeMap<bool> NodeToBoolMap;
    typedef typename Graph::template NodeMap<int> NodeToIntMap;
    typedef typename std::map<unsigned short, NodeVector> ShortToNodeVectorMap;
    typedef typename lemon::FilterNodes<const Graph, const NodeToBoolMap> FilterNodes;
    typedef typename FilterNodes::NodeIt FilteredNodeIt;
    typedef typename FilterNodes::EdgeIt FilteredEdgeIt;
    typedef typename FilterNodes::IncEdgeIt FilteredIncEdgeIt;
    typedef typename FilterNodes::NodeMap<bool> FilteredNodeToBoolMap;


    ShortVector _colors;
    LongVector _canonization;
    NodeVector _node_order;

    void init(const Molecule& mol);

    void init(const Molecule& mol, const NodeToBoolMap& filter, const Node& root);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             NodeToBoolMap& visited, ShortSet& colorSet, ShortToNodeVectorMap& colorMap,
             bool& is_tree);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             const FilterNodes& subgraph, NodeToBoolMap& visited, ShortSet& colorSet,
             ShortToNodeVectorMap& colorMap, bool& is_tree, unsigned int& node_count);

    void canonTree(const Molecule& mol, const ShortToNodeVectorMap& colorMap);

    void canonTree(const FilterNodes& subgraph, const ShortToNodeVectorMap& colorMap);

    void canonNauty(const Molecule& mol,
                    const ShortSet &colorSet,
                    const ShortToNodeVectorMap &colorMap,
                    const unsigned int atom_count);

    void canonNauty(const Molecule& mol,
                    const FilterNodes& subgraph,
                    const ShortSet &colorSet,
                    const ShortToNodeVectorMap &colorMap,
                    const Node& root,
                    const unsigned int atom_count);

  };

}

#endif //MOGLI_CANONIZATION_H
