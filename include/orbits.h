////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
//    Copyright (C) 2016-2019  Martin S. Engler                                                                       //
//                                                                                                                    //
//    This program is free software: you can redistribute it and/or modify                                            //
//    it under the terms of the GNU Lesser General Public License as published                                        //
//    by the Free Software Foundation, either version 3 of the License, or                                            //
//    (at your option) any later version.                                                                             //
//                                                                                                                    //
//    This program is distributed in the hope that it will be useful,                                                 //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                  //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                    //
//    GNU General Public License for more details.                                                                    //
//                                                                                                                    //
//    You should have received a copy of the GNU Lesser General Public License                                        //
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
