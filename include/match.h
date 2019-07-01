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

    Match(const IntToIntMap &frag_to_mol) :
        _frag_to_mol(frag_to_mol),
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

    void merge(Match &other, const IntToIntMap &isomorphism_map) {
      if (isomorphism_map.size() > 0) {
        const IntToIntMap &other_ftm = other.get_frag_to_mol();
        const IntToIntMapVector &other_mftm = other.get_merged_frag_to_mol();
        _merge_same(other_ftm, isomorphism_map);
        for (auto & it : other_mftm) {
          _merge_same(it, isomorphism_map);
        }
      }
    }

    void map_ids(const Canonization &this_canon, const Canonization &other_canon) {
      const ShortVector &this_nodes = this_canon.get_node_order();
      const ShortVector &other_nodes = other_canon.get_node_order();
      assert(this_nodes.size() == other_nodes.size());
      map_other(_frag_to_mol, this_nodes, other_nodes);
      assert(_frag_to_mol.size() == this_nodes.size());
      for (auto & it : _merged_frag_to_mol) {
        map_other(it, this_nodes, other_nodes);
      }
    }

  private:

    void _merge_same(const IntToIntMap &other, const IntToIntMap & iso_map) {
      IntToIntMap copy;
      for (auto & it : other) {
        if (iso_map.count(it.first) > 0) {
          copy[iso_map.at(it.first)] = it.second;
        }
      }
      _merged_frag_to_mol.push_back(copy);
    }

    void map_other(IntToIntMap &map, const ShortVector &this_nodes, const ShortVector &other_nodes) {
      IntToIntMap copy;
      for (int i = 0; i < this_nodes.size(); ++i) {
        if (map.count(this_nodes[i]) > 0) {
          copy[other_nodes[i]] = map.at(this_nodes[i]);
        }
      }
      map.clear();
      for (auto & it : copy) {
        map[it.first] = it.second;
      }
    }

  };

}

#endif //MOGLI_MATCH_H
