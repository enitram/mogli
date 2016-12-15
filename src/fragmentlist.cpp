//
// Created by M. Engler on 08/12/16.
//

#include "fragmentlist.h"

using namespace mogli;

bool less_vectors(const NodeVector& a,const NodeVector& b) {
  return a.size() < b.size();
}

void FragmentList::push_back_all(const Product &product, const NodeVectorVector cliques) {
  _fragments.reserve(cliques.size());
  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
    _fragments.push_back(new Fragment(product, *it));
  }
}

void FragmentList::push_back_unique(const Product &product, const NodeVectorVector cliques) {
  _fragments.reserve(cliques.size());
  std::sort(cliques.begin(), cliques.end(), less_vectors);
  std::deque<std::pair<Fragment*, Tgraph*> > deque;
  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
    Fragment* fragment = new Fragment(product, *it);
    Tgraph* tgraph = translate_graph(*fragment);
    deque.push_back(std::make_pair(fragment,tgraph));
  }
  while (deque.size() > 0) {
    std::pair<Fragment*, Tgraph*>* current = &deque.front();
    int map[current->first->get_atom_count()];
    bool iso = false;
    for (std::deque<std::pair<Fragment*, Tgraph*> >::reverse_iterator it = deque.rbegin(), end = deque.rend(); it < end-1; ++it) {
      assert(current->first->get_atom_count() <= it->first->get_atom_count());
      iso = are_subgraph_isomorphic(current->second, it->second, map);
      if (iso) {
        it->first->merge(*(current->first));
        break;
      }
    }
    if (!iso) {
      _fragments.push_back(current->first);
    }
    deque.pop_front();
  }
}
