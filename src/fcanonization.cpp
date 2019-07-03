//
// Created by martin on 7/3/19.
//

#include "fcanonization.h"

const bool mogli::FragmentCanonization::is_isomorphic(FragmentCanonization &other) const  {
  if (Canonization::is_isomorphic(other)) {
    const BoolVector & cores2 = other.get_core_nodes();

    if (_core_nodes.size() != cores2.size())
      return false;

    for (auto i1 = _core_nodes.begin(), i2 = cores2.begin(),
             ie1 = _core_nodes.end(), ie2 = cores2.end(); i1 != ie1 && i2 != ie2; ++i1, ++i2) {
      if (*i1 != *i2)
        return false;
    }

    return true;
  } else {
    return false;
  }
}