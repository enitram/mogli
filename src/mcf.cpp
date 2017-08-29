//
// Created by M. Engler on 10/01/17.
//

#include <sublad.h>
#include <numeric>
#include "../include/mcf.h"

bool less(const std::pair<int, int>& a, const std::pair<int, int>& b) {
  return a.second < b.second;
}

using namespace mogli;

void append_fragments(Product &product, NodeVectorVector& cliques, FragmentVector &fragments,
                      MatchVector &matches_mol1, MatchVector &matches_mol2) {
  for (NodeVectorVector::const_iterator it = cliques.begin(), end = cliques.end(); it != end; ++it) {
    IntToIntMap g_to_mol1, g_to_mol2;
    boost::shared_ptr<Fragment> fragment = boost::make_shared<Fragment>(product, *it, g_to_mol1, g_to_mol2);

    Match match1(g_to_mol1);
    Match match2(g_to_mol2);

    fragments.push_back(fragment);
    matches_mol1.push_back(match1);
    matches_mol2.push_back(match2);
  }
}

void run_mcf(Product &product, FragmentVector &fragments,
             MatchVector &matches_mol1, MatchVector &matches_mol2,
             int min_core_size, int max_core_size, bool maximum) {

  BronKerbosch bk(product, min_core_size, max_core_size, maximum);
  bk.run();

  NodeVectorVector cliques = bk.getMaxCliques();

  if (!maximum) {
    append_fragments(product, cliques, fragments, matches_mol1, matches_mol2);
  } else {
    int current_max = fragments.size() > 0 ? fragments[0]->get_core_atom_count() : 0;
    int cliques_size = cliques.size() > 0 ? product.get_clique_size(cliques[0]) : 0;

    if (cliques_size > current_max) {
      fragments.clear();
      matches_mol1.clear();
      matches_mol2.clear();
    }
    if (cliques_size >= current_max) {
      append_fragments(product, cliques, fragments, matches_mol1, matches_mol2);
    }
  }
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

  if (!reduce_subgraphs) {
    if (product.get_components() == 1) {
      run_mcf(product, fragments, matches_mol1, matches_mol2,
              min_core_size, max_core_size, maximum);
    } else if (!maximum) {
      for (int c = 0; c < product.get_components(); ++c) {
        Product component(product, c);
        run_mcf(component, fragments, matches_mol1, matches_mol2,
                min_core_size, max_core_size, maximum);
      }
    } else {
      // sort components descending by size
      std::vector<int> idx(product.get_components());
      std::iota(idx.begin(), idx.end(), 0);
      sort(idx.begin(), idx.end(), [&product](int i1, int i2) {return product.get_component_size(i1)
                                                                            > product.get_component_size(i2);});

      for (int c : idx) {
        // break if component smaller than current max fragment
        int current_max = fragments.size() > 0 ? fragments[0]->get_core_atom_count() : 0;
        if (product.get_component_size(c) < current_max)
          break;
        Product component(product, c);
        run_mcf(component, fragments, matches_mol1, matches_mol2,
                min_core_size, max_core_size, maximum);
      }
    }

  } else {

    std::deque<std::pair<int, int> > deque;
    FragmentVector frags;
    MatchVector matches1;
    MatchVector matches2;

    if (product.get_components() == 1) {
      run_mcf(product, frags, matches1, matches2,
              min_core_size, max_core_size, maximum);
    } else if (!maximum) {
      for (int c = 0; c < product.get_components(); ++c) {
        Product component(product, c);
        run_mcf(component, frags, matches1, matches2,
                min_core_size, max_core_size, maximum);
      }
    } else {
      // sort components descending by size
      std::vector<int> idx(product.get_components());
      std::iota(idx.begin(), idx.end(), 0);
      sort(idx.begin(), idx.end(), [&product](int i1, int i2) {return product.get_component_size(i1)
                                                                      > product.get_component_size(i2);});

      for (int c : idx) {
        // break if component smaller than current max fragment
        int current_max = frags.size() > 0 ? frags[0]->get_core_atom_count() : 0;
        if (product.get_component_size(c) < current_max)
          break;
        Product component(product, c);
        run_mcf(component, frags, matches1, matches2,
                min_core_size, max_core_size, maximum);
      }
    }

    int k = 0;
    for (auto fragment : frags) {
      deque.push_back(std::make_pair(k, (int) fragment->get_atom_count()));
      ++k;
    }

    assert(k == deque.size());

    std::sort(deque.begin(), deque.end(), less);

    fragments.reserve(deque.size());
    matches_mol1.reserve(deque.size());
    matches_mol2.reserve(deque.size());

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

void mogli::atomic_fragments(Molecule &mol, FragmentVector &fragments, MatchVector &matches, int shell) {
  typedef typename Graph::template NodeMap<NodeVector> NodeToNodeVectorMap;
  typedef std::deque<Node> NodeDeque;

  const Graph& g = mol.get_graph();
  lemon::ArcLookUp<Graph> arcLookUp(mol.get_graph());

  for (NodeIt v(g); v != lemon::INVALID; ++v) {
    NodeToBoolMap visited(mol.get_graph(), false);
    NodeToIntMap depth(mol.get_graph(), 0);
    NodeDeque queue;
    boost::shared_ptr<Fragment> fragment = boost::make_shared<Fragment>();
    Match match;

    queue.push_back(v);
    visited[v] = true;
    Node root;
    while (queue.size() > 0) {
      Node &current = queue.front();

      Node copy = fragment->add_atom(mol.get_color(current));
      match.add_frag_to_mol(fragment->get_id(copy), mol.get_id(current));
      if (current == v)
        root = copy;

      if (depth[current] < shell) {
        for (IncEdgeIt e = mol.get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
          Node w = mol.get_opposite_node(current, e);
          if (!visited[w]) {
            visited[w] = true;
            depth[w] = depth[current] + 1;
            queue.push_back(w);
          }
        }
      }
      queue.pop_front();
    }

    const Graph& gf = fragment->get_graph();
    for (NodeIt u(gf); u != lemon::INVALID; ++u) {
      for (NodeIt w = u; w != lemon::INVALID; ++w) {
        if (u == w)
          continue;

        if (arcLookUp(mol.get_node_by_id(match.frag_to_mol(fragment->get_id(u))),
                      mol.get_node_by_id(match.frag_to_mol(fragment->get_id(w)))) != lemon::INVALID) {
          fragment->add_edge(u, w);
        }
      }
    }

    fragment->set_core(root, true);
    fragment->set_shell_size(shell);

    fragments.push_back(fragment);
    matches.push_back(match);
  }

}

