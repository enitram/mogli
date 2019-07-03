#include <utility>

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

    explicit FragmentCanonization(const Fragment &fragment) : Canonization(fragment) {
      for (auto & it : _node_order) {
        _core_nodes.push_back(fragment.is_core(fragment.get_node_by_id(it)));
      }
    }

    FragmentCanonization(const ShortVector &_colors, const LongVector &_canonization,
                         const ShortVector &_node_order, BoolVector _core_nodes) :
        Canonization(_colors, _canonization, _node_order),
        _core_nodes(std::move(_core_nodes)) {}

    const BoolVector &get_core_nodes() const {
      return _core_nodes;
    }

    const bool is_isomorphic(FragmentCanonization & other) const;

  };

}


#endif //MOGLI_FCANONIZATION_H
