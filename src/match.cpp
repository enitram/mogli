//
// Created by martin on 7/3/19.
//

#include "match.h"

// public functions

void mogli::Match::merge(const Match &other, const IntToIntMap &isomorphism_map)  {
  if (!isomorphism_map.empty()) {
    const IntToIntMap &other_ftm = other.get_frag_to_mol();
    const IntToIntMapVector &other_mftm = other.get_merged_frag_to_mol();
    _merge_same(other_ftm, isomorphism_map);
    for (auto & it : other_mftm) {
      _merge_same(it, isomorphism_map);
    }
  }
}

void mogli::Match::map_ids(const Canonization &this_canon, const Canonization &other_canon)  {
  const ShortVector &this_nodes = this_canon.get_node_order();
  const ShortVector &other_nodes = other_canon.get_node_order();
  assert(this_nodes.size() == other_nodes.size());
  map_other(_frag_to_mol, this_nodes, other_nodes);
  assert(_frag_to_mol.size() == this_nodes.size());
  for (auto & it : _merged_frag_to_mol) {
    map_other(it, this_nodes, other_nodes);
  }
}

// private functions

void mogli::Match::_merge_same(const mogli::IntToIntMap &other, const mogli::IntToIntMap &iso_map)  {
  IntToIntMap copy;
  for (auto & it : other) {
    if (iso_map.count(it.first) > 0) {
      copy[iso_map.at(it.first)] = it.second;
    }
  }
  _merged_frag_to_mol.push_back(copy);
}

void mogli::Match::map_other(mogli::IntToIntMap &map, const mogli::ShortVector &this_nodes,
                             const mogli::ShortVector &other_nodes)  {
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
