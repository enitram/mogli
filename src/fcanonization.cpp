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