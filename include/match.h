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

  /**
   * Maps fragment atom IDs to molecule atom IDs.
   */
  class Match {

  private:

    IntToIntMap _frag_to_mol;
    IntToIntMapVector _merged_frag_to_mol;

  public:

    /**
     * Empty constructor.
     */
    Match() :
        _frag_to_mol(),
        _merged_frag_to_mol() {}

    /**
     * Move constructor.
     *
     * @param[in] frag_to_mol   Fragment to molecule mapping.
     */
    explicit Match(IntToIntMap frag_to_mol) :
        _frag_to_mol(std::move(frag_to_mol)),
        _merged_frag_to_mol() {
    }

    /**
     * Map fragment atom ID to molecule atom ID.
     *
     * @param[in] id    Fragment atom ID.
     * @return          Molecule atom ID.
     */
    const int frag_to_mol(const int id) const {
      if (_frag_to_mol.count(id) > 0) {
        return _frag_to_mol.at(id);
      } else {
        return -1;
      }
    }

    /**
     * Map fragment atom ID to molecule atom IDs.
     *
     * @param[in]  id   Fragment atom ID.
     * @param[out] ids  Molecule atom IDs.
     */
    void merged_frag_to_mol(const int id, IntVector &ids) const {
      for (auto & it : _merged_frag_to_mol) {
        if (it.count(id) > 0) {
          ids.push_back(it.at(id));
        }
      }
    }

    /**
     * Returns the fragment to molecule atom mapping.
     *
     * @return  Fragment to molecule atom mapping.
     */
    const IntToIntMap& get_frag_to_mol() const {
      return _frag_to_mol;
    }

    /**
     * Returns all mapped molecule IDs.
     *
     * @param[out] ids  Mapped molecule IDs.
     */
    void get_atom_ids(IntVector& ids) const {
      for (auto el : _frag_to_mol) {
        ids.push_back(el.second);
      }
    }

    /**
     * Returns the fragment to molecule atom mappings.
     *
     * @return  Fragment to molecule atom mappings.
     */
    const IntToIntMapVector& get_merged_frag_to_mol() const {
      return _merged_frag_to_mol;
    }

    /**
     * Add a new fragment to molecule atom mapping.
     *
     * @param[in] from  Fragment atom ID.
     * @param[in] to    Molecule atom ID.
     */
    void add_frag_to_mol(int from, int to) {
      _frag_to_mol[from] = to;
    }

    /**
     * Add new fragment to molecule atom mappings.
     *
     * @param[in] ftm   Fragment to molecule atom mappings.
     */
    void add_merged_frag_to_mol(IntToIntMap &ftm) {
      _merged_frag_to_mol.push_back(ftm);
    }

    /**
     * Merge with another match object.
     *
     * @param[in] other             Other match object.
     * @param[in] isomorphism_map   Isomorphism map.
     */
    void merge(const Match &other, const IntToIntMap &isomorphism_map);

    /**
     * Transform fragment to molecule atom mapping to match another molecular graph.
     *
     * @param[in] this_canon    Canonization of the molecular graph this match object maps to.
     * @param[in] other_canon   Canonization of another molecular graph.
     */
    void map_ids(const Canonization &this_canon, const Canonization &other_canon);

  private:

    void _merge_same(const IntToIntMap &other, const IntToIntMap & iso_map);

    void map_other(IntToIntMap &map, const ShortVector &this_nodes, const ShortVector &other_nodes);

  };

}

#endif //MOGLI_MATCH_H
