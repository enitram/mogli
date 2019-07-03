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

#ifndef MOGLI_MCF_H
#define MOGLI_MCF_H

#include "match.h"
#include "bronkerbosch.h"
#include "subgraph_isomorphism.h"

namespace mogli {

  typedef SharedPtrVector<Fragment>::type FragmentVector;
  typedef std::vector<Match> MatchVector;

  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, unsigned int min_core_size, unsigned int max_core_size,
                                Product::GenerationType prod_gen,
                                bool reduce_subgraphs, bool maximum,
                                int timeout_seconds);

  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, unsigned int min_core_size,
                                Product::GenerationType prod_gen,
                                bool reduce_subgraphs, bool maximum,
                                int timeout_seconds);

  void atomic_fragments(Molecule &mol, FragmentVector &fragments, MatchVector &matches, int shell);

}

#endif //MOGLI_MCF_H