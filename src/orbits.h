//
// Created by M. Engler on 02/12/16.
//

#ifndef MOGLI_ORBITS_H
#define MOGLI_ORBITS_H


#include "molecule.h"
#include <nauty.h>
#include <malloc.h>

namespace mogli {

  class Orbits {

  public:

    Orbits(const Molecule &mol) : _orbits(mol.get_graph(), -1) {
      init(mol);
    }

    int get_orbit_id(const Node &v) {
      return _orbits[v];
    }

    bool same_orbit(const Node &u, const Node &v) {
      return _orbits[u] == _orbits[v];
    }

  private:

    typedef std::set<unsigned short> ShortSet;
    typedef typename std::map<unsigned short, NodeVector> ShortToNodeVectorMap;
    typedef typename Graph::template NodeMap<int> NodeToIntMap;

    NodeToIntMap _orbits;

    void init(const Molecule &mol);

    void dfs(const Node& current,
             const Molecule& mol,
             NodeToBoolMap& visited,
             ShortSet& colorSet,
             ShortToNodeVectorMap& colorMap);

    void orbitsNauty(const Molecule& mol,
                    const ShortSet &colorSet,
                    const ShortToNodeVectorMap &colorMap,
                    const unsigned int atom_count);

  };

}




#endif //MOGLI_ORBITS_H
