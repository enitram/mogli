#include <utility>

//
// Created by M. Engler on 19/01/17.
//

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
