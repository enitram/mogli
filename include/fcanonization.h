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
