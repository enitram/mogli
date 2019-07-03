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

#ifndef MOGLI_SUBGRAPH_ISOMORPHISM_H
#define MOGLI_SUBGRAPH_ISOMORPHISM_H

#include <sublad.h>
#include "molecule.h"


namespace mogli {

  Tgraph* translate_graph(const Molecule &mol, IntVector &node_ids);
  
  void free_graph(Tgraph* graph);

  void translate_maps(const IntVector &node_ids_small, const IntVector &node_ids_large,
                      const int in_iso_map[], IntToIntMap &out_iso_map);

  bool are_subgraph_isomorphic(const Molecule &mol_small, const Molecule &mol_large,
                               IntToIntMap isomorphism_map);

  bool are_subgraph_isomorphic(Tgraph* graph_small, Tgraph* graph_large, int isomorphism_map[]);

}

#endif //MOGLI_SUBGRAPH_ISOMORPHISM_H
