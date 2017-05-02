//
// Created by M. Engler on 10/01/17.
//

#include <sublad.h>
#include "mcf.h"

bool less(const std::pair<int, int>& a, const std::pair<int, int>& b) {
  return a.second < b.second;
}

void mogli::maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                     FragmentVector &fragments,
                                     MatchVector &matches_mol1, MatchVector &matches_mol2,
                                     int shell, unsigned int min_core_size,
                                     Product::GenerationType prod_gen, bool reduce_subgraphs, bool maximum) {
  maximal_common_fragments(mol1, mol2, fragments, matches_mol1, matches_mol2, shell, min_core_size,
                           std::numeric_limits<int>::max(), prod_gen, reduce_subgraphs, maximum);
}

void mogli::maximal_common_fragments(Molecule &mol1, Molecule &mol2,
                                     FragmentVector &fragments,
                                     MatchVector &matches_mol1, MatchVector &matches_mol2,
                                     int shell, unsigned int min_core_size, unsigned int max_core_size,
                                     Product::GenerationType prod_gen, bool reduce_subgraphs, bool maximum) {

  Product product(mol1, mol2, shell, prod_gen, min_core_size, max_core_size);

//  std::ofstream ofs("/home/martin/workspace/mogli/product.dot", std::ifstream::out);
//  product.print_dot(ofs);
//  ofs.close();

  BronKerbosch bk(product, min_core_size, max_core_size, maximum);
  bk.run();

  NodeVectorVector cliques = bk.getMaxCliques();
  if (!reduce_subgraphs) {

    for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
      IntToIntMap g_to_mol1, g_to_mol2;
      boost::shared_ptr<Fragment> fragment = boost::make_shared<Fragment>(product, *it, g_to_mol1, g_to_mol2);

      if (fragment->get_core_atom_count() >= min_core_size && fragment->get_core_atom_count() <= max_core_size) {
        Match match1(g_to_mol1);
        Match match2(g_to_mol2);

        fragments.push_back(fragment);
        matches_mol1.push_back(match1);
        matches_mol2.push_back(match2);
      }
    }

  } else {

    std::deque<std::pair<int, int> > deque;
    FragmentVector frags;
    MatchVector matches1;
    MatchVector matches2;

    int k = 0;
    for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
      IntToIntMap g_to_mol1, g_to_mol2;
      boost::shared_ptr<Fragment> fragment = boost::make_shared<Fragment>(product, *it, g_to_mol1, g_to_mol2);

      if (fragment->get_core_atom_count() >= min_core_size && fragment->get_core_atom_count() <= max_core_size) {
        Match match1(g_to_mol1);
        Match match2(g_to_mol2);

        frags.push_back(fragment);
        matches1.push_back(match1);
        matches2.push_back(match2);
        deque.push_back(std::make_pair(k, fragment->get_atom_count()));
        ++k;
      }
    }

    assert(k == deque.size());

    std::sort(deque.begin(), deque.end(), less);

    fragments.reserve(cliques.size());
    matches_mol1.reserve(cliques.size());
    matches_mol2.reserve(cliques.size());

    while (deque.size() > 0) {

      int current = deque.begin()->first;

      assert(0 <= current);
      assert(current <= frags.size());
      assert(current <= matches1.size());
      assert(current <= matches2.size());

      IntVector node_ids_current;
      node_ids_current.reserve(frags.at(current)->get_atom_count());
      Tgraph* graph_small = translate_graph(*frags.at(current), node_ids_current);
      int n = frags.at(current)->get_atom_count();
      int map[n];
      bool are_sub_iso = false;

      for (std::deque<std::pair<int, int> >::reverse_iterator it = deque.rbegin(), end = deque.rend(); it < end-1; ++it) {

        int other = it->first;

        assert(0 <= other);
        assert(other <= frags.size());
        assert(other <= matches1.size());
        assert(other <= matches2.size());

        assert(current != other);
        assert(current == deque.begin()->first);
        assert(other == it->first);
        assert(frags.at(current)->get_atom_count() <= frags.at(other)->get_atom_count());

        IntVector node_ids_other;
        node_ids_other.reserve(frags.at(other)->get_atom_count());
        Tgraph* graph_large = translate_graph(*frags.at(other), node_ids_other);
        std::fill_n(map, n, -1);

        are_sub_iso = are_subgraph_isomorphic(graph_small, graph_large, map);
        free_graph(graph_large);

        if (are_sub_iso) {
          IntToIntMap isomorphism_map;
          translate_maps(node_ids_current, node_ids_other, map, isomorphism_map);
          matches1.at(other).merge(matches1.at(current), isomorphism_map);
          matches2.at(other).merge(matches2.at(current), isomorphism_map);
          break;
        }

      }

      if (!are_sub_iso) {
        fragments.push_back(frags.at(current));
        matches_mol1.push_back(matches1.at(current));
        matches_mol2.push_back(matches2.at(current));
      }

      free_graph(graph_small);
      deque.pop_front();

    }
  }

}
