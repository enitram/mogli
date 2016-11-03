//
// Created by M. Engler on 10/20/16.
//

#ifndef MOGLI_CANONIZATION_H
#define MOGLI_CANONIZATION_H

#include "molecule.h"
#include <nauty.h>
#include <malloc.h>
#include <bitset>

namespace mogli {

  typedef std::vector<unsigned short> ShortVector;
  typedef std::vector<unsigned long> LongVector;

  class Canonization {

  public:

    Canonization() {}

    Canonization(const Molecule& mol) : _colors(), _canonization() {
      init(mol);
    }

    Canonization(const Molecule& mol, const NodeToBoolMap& filter) {
      init(mol, filter);
    }

    Canonization(const ShortVector &_colors, const LongVector &_canonization) : _colors(_colors),
                                                                             _canonization(_canonization) {}
    const ShortVector &get_colors() const {
      return _colors;
    }

    const LongVector &get_canonization() const {
      return _canonization;
    }

  private:

    typedef std::set<unsigned short> ShortSet;
    typedef typename Graph::template NodeMap<bool> NodeToBoolMap;
    typedef typename std::vector<Node> NodeVector;
    typedef typename std::map<unsigned short, NodeVector> ShortToNodeVectorMap;
    typedef typename lemon::FilterNodes<const Graph, const NodeToBoolMap> FilterNodes;
    typedef typename FilterNodes::NodeIt FilteredNodeIt;
    typedef typename FilterNodes::EdgeIt FilteredEdgeIt;
    typedef typename FilterNodes::IncEdgeIt FilteredIncEdgeIt;
    typedef typename FilterNodes::NodeMap<bool> FilteredNodeToBoolMap;

    ShortVector _colors;
    LongVector _canonization;

    void init(const Molecule& mol);

    void init(const Molecule& mol, const NodeToBoolMap& filter);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             NodeToBoolMap& visited, ShortSet& colorSet, ShortToNodeVectorMap& colorMap,
             bool& is_tree);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             const FilterNodes& subgraph, NodeToBoolMap& visited, ShortSet& colorSet,
             ShortToNodeVectorMap& colorMap, bool& is_tree);

    void canonTree(const Molecule& mol, const ShortToNodeVectorMap& colorMap);

    void canonTree(const FilterNodes& subgraph, const ShortToNodeVectorMap& colorMap);

    void canonNauty(const Molecule& mol, const ShortToNodeVectorMap &colorMap, const unsigned int atom_count);

    void canonNauty(const FilterNodes& subgraph, const ShortToNodeVectorMap &colorMap, const unsigned int atom_count);

  };

}

#endif //MOGLI_CANONIZATION_H
