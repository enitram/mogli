//
// Created by M. Engler on 10/01/17.
//

# include "mcf.h"

bool mogli::less(const std::pair<Fragment*, Tgraph*>& a,const std::pair<Fragment*, Tgraph*>& b) {
  return a.first->get_atom_count() < b.first->get_atom_count();
}

void mogli::maximal_common_fragments(Molecule &mol1, Molecule &mol2, FragmentVector &fragments,
                                     int shell, Product::GenerationType prod_gen, std::string unique_node_property) {

  Product product(mol1, mol2, shell, prod_gen);

  BronKerbosch bk(product);
  bk.run();

  NodeVectorVector cliques = NodeVectorVector(bk.getMaxCliques());

  std::deque<std::pair<Fragment*, Tgraph*> > deque;
  int i = 0;
  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
    Fragment* fragment = new Fragment(product, *it, unique_node_property);
    Tgraph* tgraph = translate_graph(*fragment);
    deque.push_back(std::make_pair(fragment,tgraph));
  }

  std::sort(deque.begin(), deque.end(), less);

  fragments.reserve(cliques.size());
  while (deque.size() > 0) {
    std::pair<Fragment*, Tgraph*>* current = &deque.front();
    int map[current->first->get_atom_count()];
    bool iso = false;
    for (std::deque<std::pair<Fragment*, Tgraph*> >::reverse_iterator it = deque.rbegin(), end = deque.rend(); it < end-1; ++it) {
      assert(current->first->get_atom_count() <= it->first->get_atom_count());
      iso = are_subgraph_isomorphic(current->second, it->second, map);
      if (iso) {
        it->first->merge(*(current->first));
        delete current->first;
        break;
      }
    }
    if (!iso) {
      fragments.push_back(current->first);
    }
    free_graph(current->second);
    deque.pop_front();
  }

}
