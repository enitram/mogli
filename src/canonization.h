//
// Created by martin on 10/20/16.
//

#ifndef MOGLI_CANONIZATION_H
#define MOGLI_CANONIZATION_H

#include "molecule.h"
#include <nauty.h>
#include <malloc.h>
#include <bitset>

namespace mogli {

  typedef std::set<unsigned short> ShortSet;
  typedef std::vector<unsigned long> LongVector;

  class Canonization {

  public:

    Canonization(Molecule& mol) : _colors(), _canonization() {
      init(mol);
    }

    virtual ~Canonization() {}

    const ShortSet &get_colors() const {
      return _colors;
    }

    const LongVector &get_canonization() const {
      return _canonization;
    }

  private:

    typedef typename Graph::template NodeMap<bool> NodeToBoolMap;
    typedef typename std::vector<Node> NodeVector;
    typedef typename std::map<unsigned short, NodeVector> ShortToNodeVectorMap;

    ShortSet _colors;
    LongVector _canonization;

    void init(const Molecule& mol);

    void dfs(const Node& current, const Node& last,
             const Molecule& mol, const Graph& g, NodeToBoolMap& visited,
             ShortToNodeVectorMap& colorMap, bool& is_tree);

    void canonTree(const Graph& g, const ShortToNodeVectorMap& colorMap);

    void canonNauty(const Graph& g, const ShortToNodeVectorMap &colorMap, const unsigned int atom_count);

  };

}

#endif //MOGLI_CANONIZATION_H
