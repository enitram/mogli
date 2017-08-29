//
// Created by M. Engler on 10/01/17.
//

#ifndef MOGLI_MCF_H
#define MOGLI_MCF_H

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "match.h"
#include "bronkerbosch.h"
#include "subgraph_isomorphism.h"

namespace mogli {

  typedef std::vector<boost::shared_ptr<Fragment> > FragmentVector;
  typedef std::vector<Match> MatchVector;

  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, unsigned int min_core_size, unsigned int max_core_size,
                                Product::GenerationType prod_gen,
                                bool reduce_subgraphs, bool maximum);

  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, unsigned int min_core_size,
                                Product::GenerationType prod_gen,
                                bool reduce_subgraphs, bool maximum);

  void atomic_fragments(Molecule &mol, FragmentVector &fragments, MatchVector &matches, int shell);

}

#endif //MOGLI_MCF_H