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

//  bool less(const std::pair<Fragment*, Tgraph*>& a, const std::pair<Fragment*, Tgraph*>& b);

  void maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                FragmentVector &fragments,
                                MatchVector &matches_mol1, MatchVector &matches_mol2,
                                int shell, int core_size_limit,
                                Product::GenerationType prod_gen, bool reduce_subgraphs);

}

#endif //MOGLI_MCF_H