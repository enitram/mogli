//
// Created by M. Engler on 10/01/17.
//

#ifndef MOGLI_MCF_H
#define MOGLI_MCF_H

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "fragment.h"
#include "bronkerbosch.h"
#include "subgraph_isomorphism.h"

namespace mogli {

  typedef boost::ptr_vector<Fragment> FragmentVector;

  bool less(const std::pair<Fragment*, Tgraph*>& a, const std::pair<Fragment*, Tgraph*>& b);

  void maximal_common_fragments(Molecule &mol1, Molecule &mol2, FragmentVector &fragments,
                                int shell, Product::GenerationType prod_gen, std::string unique_node_property);

}

#endif //MOGLI_MCF_H