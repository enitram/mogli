//
// Created by M. Engler on 10/20/16.
//

#ifndef MOGLI_CANONIZATION_H
#define MOGLI_CANONIZATION_H

#include "molecule.h"
#include <nauty.h>
#include <malloc.h>

namespace mogli {


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

    Canonization(const ShortVector &_colors, const LongVector &_canonization, const ShortVector &_node_order) :
        _colors(_colors), _canonization(_canonization), _node_order(_node_order) {}

    const ShortVector &get_colors() const {
      return _colors;
    }

    const LongVector &get_canonization() const {
      return _canonization;
    }

    const ShortVector &get_node_order() const {
      return _node_order;
    }

    const bool is_isomorphic(Canonization &other) const {
      const ShortVector& colors2 = other.get_colors();

      if (_colors.size() != colors2.size())
        return false;

      const LongVector& canonization2 = other.get_canonization();

      if (_canonization.size() != canonization2.size())
        return false;

      for (ShortVector::const_iterator i1 = _colors.begin(), i2 = colors2.begin(),
               ie1 = _colors.end(), ie2 = colors2.end();
           i1 != ie1 && i2 != ie2; ++i1, ++i2) {
        if (*i1 != *i2)
          return false;
      }

      for (LongVector::const_iterator i1 = _canonization.begin(), i2 = canonization2.begin(),
               ie1 = _canonization.end(), ie2 = canonization2.end();
           i1 != ie1 && i2 != ie2; ++i1, ++i2) {
        if (*i1 != *i2)
          return false;
      }

      return true;
    }

  protected:

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
    ShortVector _node_order;

    void init(const Molecule& mol);

    void init(const Molecule& mol, const NodeToBoolMap& filter, const Node& root);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             NodeToBoolMap& visited, ShortSet& colorSet, ShortToNodeVectorMap& colorMap);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             const FilterNodes& subgraph, NodeToBoolMap& visited, ShortSet& colorSet,
             ShortToNodeVectorMap& colorMap, unsigned int& node_count);

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
