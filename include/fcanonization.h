//
// Created by M. Engler on 30/01/17.
//

#ifndef MOGLI_FCANONIZATION_H
#define MOGLI_FCANONIZATION_H

#include "canonization.h"
#include "fragment.h"

namespace mogli {

  class FragmentCanonization : public Canonization {
  private:

    BoolVector _core_nodes;

  public:

    FragmentCanonization() : Canonization() {}

    FragmentCanonization(const Fragment &fragment) : Canonization(fragment) {
      for (auto & it : _node_order) {
        _core_nodes.push_back(fragment.is_core(fragment.get_node_by_id(it)));
      }
    }

    FragmentCanonization(const ShortVector &_colors, const LongVector &_canonization,
                         const ShortVector &_node_order, const BoolVector &_core_nodes) :
        Canonization(_colors, _canonization, _node_order),
        _core_nodes(_core_nodes) {}

    const BoolVector &get_core_nodes() const {
      return _core_nodes;
    }

    const bool is_isomorphic(FragmentCanonization & other) const {
      if (Canonization::is_isomorphic(other)) {
        BoolVector cores2 = other.get_core_nodes();

        if (_core_nodes.size() != cores2.size())
          return false;

        for (BoolVector::const_iterator i1 = _core_nodes.begin(), i2 = cores2.begin(),
                 ie1 = _core_nodes.end(), ie2 = cores2.end(); i1 != ie1 && i2 != ie2; ++i1, ++i2) {
          if (*i1 != *i2)
            return false;
        }

        return true;
      } else {
        return false;
      }
    }

  };

}


#endif //MOGLI_FCANONIZATION_H
