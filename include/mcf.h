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

  /**
   * @brief Computes maximal common fragments of two molecular graphs.
   *
   * The heart and soul of this library. See <a href="https://doi.org/10.7287/peerj.preprints.3250v1">this paper</a>
   * for more information. The data reduction rule with the most speedup is Product::GenerationType::UNCON_DEG_1. It is
   * recommended to always use this rule, the other rules are mainly for evaluation.
   *
   * @param[in]  mol1               First molecular graph.
   * @param[in]  mol2               Second molecular graph.
   * @param[out] fragments          Maximal common fragments.
   * @param[out] matches_mol1       Match objects mapping from fragments to the first molecular graph.
   * @param[out] matches_mol2       Match objects mapping from fragments to the second molecular graph.
   * @param[in]  shell              Shell size. Maximal number of bonds from any core atom in the fragments.
   * @param[in]  min_core_size      Minimal number of core atoms for each fragment.
   * @param[in]  max_core_size      Maximal number of core atoms for each fragment.
   * @param[in]  prod_gen           Product graph data reduction rule.
   * @param[in]  reduce_subgraphs   Merges resulting fragments, if one is a subgraph of the other.
   * @param[in]  maximum            If true, reports only the largest fragments.
   * @param[in]  timeout_seconds    Timeout in seconds.
   */
  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, unsigned int min_core_size, unsigned int max_core_size,
                                Product::GenerationType prod_gen,
                                bool reduce_subgraphs, bool maximum,
                                int timeout_seconds);
  /**
   * @brief Computes maximal common fragments of two molecular graphs.
   *
   * The heart and soul of this library. See <a href="https://doi.org/10.7287/peerj.preprints.3250v1">this paper</a>
   * for more information. The data reduction rule with the most speedup is Product::GenerationType::UNCON_DEG_1. It is
   * recommended to always use this rule, the other rules are mainly for evaluation.
   *
   * @param[in]  mol1               First molecular graph.
   * @param[in]  mol2               Second molecular graph.
   * @param[out] fragments          Maximal common fragments.
   * @param[out] matches_mol1       Match objects mapping from fragments to the first molecular graph.
   * @param[out] matches_mol2       Match objects mapping from fragments to the second molecular graph.
   * @param[in]  shell              Shell size. Maximal number of bonds from any core atom in the fragments.
   * @param[in]  min_core_size      Minimal number of core atoms for each fragment.
   * @param[in]  prod_gen           Product graph data reduction rule.
   * @param[in]  reduce_subgraphs   Merges resulting fragments, if one is a subgraph of the other.
   * @param[in]  maximum            If true, reports only the largest fragments.
   * @param[in]  timeout_seconds    Timeout in seconds.
   */
  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, unsigned int min_core_size,
                                Product::GenerationType prod_gen,
                                bool reduce_subgraphs, bool maximum,
                                int timeout_seconds);

  /**
   * Iterates all atoms of a molecular graph and returns them as fragments with a single-atom core.
   *
   * @param[in]  mol          Molecular graph.
   * @param[out] fragments    Atomic fragments.
   * @param[out] matches      Match objects mapping from fragments to the molecular graph.
   * @param[in]  shell        Shell size of the subgraphs. Maximal number of bonds from the center atom.
   */
  void atomic_fragments(Molecule &mol, FragmentVector &fragments, MatchVector &matches, int shell);

}

#endif //MOGLI_MCF_H