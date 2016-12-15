//
// Created by M. Engler on 08/12/16.
//

#ifndef MOGLI_FRAGMENTLIST_H
#define MOGLI_FRAGMENTLIST_H

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "fragment.h"
#include "subgraph_isomorphism.h"

namespace mogli {

  // TODO test this, make python bindings, change FDB

  typedef boost::void_ptr_iterator<std::vector<void *>::const_iterator, const Fragment> FragmentIterator;

  class FragmentList {

  private:

    typedef boost::ptr_vector<Fragment> FragmentVector;

    FragmentVector _fragments;

  public:

    FragmentList() : _fragments() {}

    FragmentIterator begin() {
      return _fragments.begin();
    }

    FragmentIterator end() {
      return _fragments.end();
    }

    void push_back_all(const Product &product, const NodeVectorVector cliques);

    void push_back_unique(const Product &product, const NodeVectorVector cliques);

  };

}


#endif //MOGLI_FRAGMENTLIST_H
