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

#ifndef MOGLI_MATCH_H
#define MOGLI_MATCH_H

#include "fragment.h"

namespace mogli {

  class Match {

  private:

    IntToIntMap _frag_to_mol;
    IntToIntMapVector _merged_frag_to_mol;

  public:

    Match() :
        _frag_to_mol(),
        _merged_frag_to_mol() {}

    explicit Match(IntToIntMap frag_to_mol) :
        _frag_to_mol(std::move(frag_to_mol)),
        _merged_frag_to_mol() {
    }

    const int frag_to_mol(const int id) const {
      if (_frag_to_mol.count(id) > 0) {
        return _frag_to_mol.at(id);
      } else {
        return -1;
      }
    }

    void merged_frag_to_mol(const int id, IntVector &ids) const {
      for (auto & it : _merged_frag_to_mol) {
        if (it.count(id) > 0) {
          ids.push_back(it.at(id));
        }
      }
    }

    const IntToIntMap& get_frag_to_mol() const {
      return _frag_to_mol;
    }

    void get_atom_ids(IntVector& ids) const {
      for (auto el : _frag_to_mol) {
        ids.push_back(el.second);
      }
    }

    const IntToIntMapVector& get_merged_frag_to_mol() const {
      return _merged_frag_to_mol;
    }

    void add_frag_to_mol(int from, int to) {
      _frag_to_mol[from] = to;
    }

    void add_merged_frag_to_mol(IntToIntMap &ftm) {
      _merged_frag_to_mol.push_back(ftm);
    }

    void merge(const Match &other, const IntToIntMap &isomorphism_map);

    void map_ids(const Canonization &this_canon, const Canonization &other_canon);

  private:

    void _merge_same(const IntToIntMap &other, const IntToIntMap & iso_map);

    void map_other(IntToIntMap &map, const ShortVector &this_nodes, const ShortVector &other_nodes);

  };

}

#endif //MOGLI_MATCH_H
