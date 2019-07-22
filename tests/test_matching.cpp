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

#include "catch2/catch.hpp"
#include "mcf.h"
#include "test_fixtures.h"

using namespace mogli;


TEST_CASE("lookup", "[algo]") {

  Graph graph;
  auto u = graph.addNode();
  auto v = graph.addNode();
  auto e = graph.addEdge(u, v);

  lemon::ArcLookUp<Graph> lookup(graph);

  REQUIRE(lookup(u, v) != lemon::INVALID);
  REQUIRE(lookup(v, u) != lemon::INVALID);

  auto arc = lookup(u, v);
  REQUIRE(graph.source(arc) == graph.u(e));
  REQUIRE(graph.target(arc) == graph.v(e));

  arc = lookup(v, u);
  REQUIRE(graph.source(arc) == graph.v(e));
  REQUIRE(graph.target(arc) == graph.u(e));

}

TEST_CASE("product_noopt", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);

  Product p0(mol1, mol2, 0, Product::GenerationType::NO_OPT, 0);
  Product p1(mol1, mol2, 1, Product::GenerationType::NO_OPT, 0);
  Product p2(mol1, mol2, 2, Product::GenerationType::NO_OPT, 0);

  const auto & g0 = p0.get_graph();
  const auto & g1 = p1.get_graph();
  const auto & g2 = p2.get_graph();

  // product nodes: 2 * 2 C-nodes + 6 * 6 H-nodes
  REQUIRE(lemon::countNodes(g0) == 2 * 2 + 6 * 6 );
  REQUIRE(lemon::countNodes(g1) == 2 * 2 + 6 * 6 );
  REQUIRE(lemon::countNodes(g2) == 2 * 2 + 6 * 6 );

  int g0_edges = 0, g0_c_edges = 0;
  for (auto e = EdgeIt(g0); e != lemon::INVALID; ++e, ++g0_edges) {
    if (p0.is_connectivity_edge(e)) {
      ++g0_c_edges;
    }
  }

  int g1_edges = 0, g1_c_edges = 0;
  for (auto e = EdgeIt(g1); e != lemon::INVALID; ++e, ++g1_edges) {
    if (p1.is_connectivity_edge(e)) {
      ++g1_c_edges;
    }
  }

  int g2_edges = 0, g2_c_edges = 0;
  for (auto e = EdgeIt(g2); e != lemon::INVALID; ++e, ++g2_edges) {
    if (p2.is_connectivity_edge(e)) {
      ++g2_c_edges;
    }
  }

  // c-edges: 2 between the C-product nodes and for each C-node: 3 * 3 possible mappings of the adjacent H-nodes,
  REQUIRE(g0_c_edges == 4 * 9 + 2);
  REQUIRE(g1_c_edges == 4 * 9 + 2);
  REQUIRE(g2_c_edges == 4 * 9 + 2);

  REQUIRE(g0_edges == 524);
  REQUIRE(g1_edges == 524);
  REQUIRE(g2_edges == 524);

}

TEST_CASE("product_deg1", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);

  Product p(mol1, mol2, 1, Product::GenerationType::DEG_1, 0);

  const auto & g = p.get_graph();

  // product nodes: 2 * 2 C-nodes
  REQUIRE(lemon::countNodes(g) == 4);

  // 2 possible mappings (c-edges)
  REQUIRE(lemon::countEdges(g) == 2);
  for (auto e = EdgeIt(g); e != lemon::INVALID; ++e) {
    REQUIRE(p.is_connectivity_edge(e));
  }

}


TEST_CASE("product_uncon", "[algo]") {

  Molecule mol1, mol2;
  LGFIOConfig config("label", "atomType");

  mol1.read_lgf(GENERIC_1, config);
  mol2.read_lgf(GENERIC_2, config);

  Product p(mol1, mol2, 1, Product::GenerationType::UNCON, 1);

  const auto & g = p.get_graph();

  REQUIRE(p.get_components() == 2);

  Product _p1(p, 0), _p2(p, 1);

  Product & p1 = (lemon::countNodes(_p1.get_graph()) < lemon::countNodes(_p2.get_graph())) ? _p1 : _p2;
  Product & p2 = (lemon::countNodes(_p1.get_graph()) < lemon::countNodes(_p2.get_graph())) ? _p2 : _p1;

  const auto & g1 = p1.get_graph();
  const auto & g2 = p2.get_graph();

  REQUIRE(lemon::countNodes(g1) == 4);
  REQUIRE(lemon::countNodes(g2) == 5);

  REQUIRE(lemon::countEdges(g1) == 6);
  REQUIRE(lemon::countEdges(g2) == 6);

  int g1_c_edges = 0, g2_c_edges = 0;
  for (auto e = EdgeIt(g1); e != lemon::INVALID; ++e) {
    if (p1.is_connectivity_edge(e)) {
      ++g1_c_edges;
    }
  }
  for (auto e = EdgeIt(g2); e != lemon::INVALID; ++e) {
    if (p2.is_connectivity_edge(e)) {
      ++g2_c_edges;
    }
  }

  REQUIRE(g1_c_edges == 3);
  REQUIRE(g2_c_edges == 4);

}

TEST_CASE("product_uncon_deg1", "[algo]") {

  Molecule mol1, mol2;
  LGFIOConfig config("label", "atomType");

  mol1.read_lgf(GENERIC_1, config);
  mol2.read_lgf(GENERIC_2, config);

  Product p(mol1, mol2, 1, Product::GenerationType::UNCON_DEG_1, 1);

  const auto & g = p.get_graph();

  REQUIRE(p.get_components() == 2);

  Product _p1(p, 0), _p2(p, 1);

  Product & p1 = (lemon::countNodes(_p1.get_graph()) < lemon::countNodes(_p2.get_graph())) ? _p1 : _p2;
  Product & p2 = (lemon::countNodes(_p1.get_graph()) < lemon::countNodes(_p2.get_graph())) ? _p2 : _p1;

  const auto & g1 = p1.get_graph();
  const auto & g2 = p2.get_graph();

  REQUIRE(lemon::countNodes(g1) == 1);
  REQUIRE(lemon::countNodes(g2) == 2);

  REQUIRE(lemon::countEdges(g1) == 0);
  REQUIRE(lemon::countEdges(g2) == 1);

  for (auto e = EdgeIt(g2); e != lemon::INVALID; ++e) {
    REQUIRE(p2.is_connectivity_edge(e));
  }

}


TEST_CASE("bronkerbosch", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);

  Product p(mol1, mol2, 1, Product::GenerationType::DEG_1, 0);

  BronKerbosch bk(p, 0, std::numeric_limits<int>::max(), false);
  bk.run(TIMEOUT);

  auto cliques = bk.getMaxCliques();

  REQUIRE(cliques.size() == 2);
  REQUIRE(cliques[0].size() == 2);
  REQUIRE(cliques[1].size() == 2);

}

TEST_CASE("mcf_isomorphic_graphs", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);

  FragmentVector frag_noopt, frag_deg1, frag_uncon, frag_uncon_deg1;
  MatchVector matches1, matches2;

  maximal_common_fragments(
      mol1, mol2, frag_noopt, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);
  maximal_common_fragments(
      mol1, mol2, frag_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::DEG_1,
      false, false, TIMEOUT);
  maximal_common_fragments(
      mol1, mol2, frag_uncon, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON,
      false, false, TIMEOUT);
  maximal_common_fragments(
      mol1, mol2, frag_uncon_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON_DEG_1,
      false, false, TIMEOUT);

  REQUIRE(!frag_noopt.empty());
  REQUIRE(!frag_deg1.empty());
  REQUIRE(!frag_uncon.empty());
  REQUIRE(!frag_uncon_deg1.empty());

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag_noopt.begin(), frag_noopt.end(), comparator);
  std::sort(frag_deg1.begin(), frag_deg1.end(), comparator);
  std::sort(frag_uncon.begin(), frag_uncon.end(), comparator);
  std::sort(frag_uncon_deg1.begin(), frag_uncon_deg1.end(), comparator);

  REQUIRE(frag_noopt[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_deg1[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon_deg1[0]->get_atom_count() == mol1.get_atom_count());

}

TEST_CASE("mcf_isomorphic_graphs_big", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(MOLID3246);
  mol2.read_lgf(MOLID3246);

  FragmentVector frag_noopt, frag_deg1, frag_uncon, frag_uncon_deg1;
  MatchVector matches1, matches2;

  maximal_common_fragments(
      mol1, mol2, frag_noopt, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::DEG_1,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_uncon, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_uncon_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON_DEG_1,
      false, false, TIMEOUT_BIG);

  REQUIRE(!frag_noopt.empty());
  REQUIRE(!frag_deg1.empty());
  REQUIRE(!frag_uncon.empty());
  REQUIRE(!frag_uncon_deg1.empty());

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag_noopt.begin(), frag_noopt.end(), comparator);
  std::sort(frag_deg1.begin(), frag_deg1.end(), comparator);
  std::sort(frag_uncon.begin(), frag_uncon.end(), comparator);
  std::sort(frag_uncon_deg1.begin(), frag_uncon_deg1.end(), comparator);

  REQUIRE(frag_noopt[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_deg1[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon_deg1[0]->get_atom_count() == mol1.get_atom_count());

}

TEST_CASE("mcf_isomorphic_graphs_large", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(PACLITAXEL);
  mol2.read_lgf(PACLITAXEL);

  FragmentVector frag_noopt, frag_deg1, frag_uncon, frag_uncon_deg1;
  MatchVector matches1, matches2;

  maximal_common_fragments(
      mol1, mol2, frag_uncon, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON,
      false, false, TIMEOUT_LARGE);
  maximal_common_fragments(
      mol1, mol2, frag_uncon_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON_DEG_1,
      false, false, TIMEOUT_LARGE);

  REQUIRE(!frag_uncon.empty());
  REQUIRE(!frag_uncon_deg1.empty());

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag_uncon.begin(), frag_uncon.end(), comparator);
  std::sort(frag_uncon_deg1.begin(), frag_uncon_deg1.end(), comparator);

  REQUIRE(frag_uncon[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon_deg1[0]->get_atom_count() == mol1.get_atom_count());

}

TEST_CASE("mcf_subisomorphic_graphs", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHYL);
  mol2.read_lgf(ETHANE_2);

  FragmentVector frag_noopt, frag_deg1, frag_uncon, frag_uncon_deg1;
  MatchVector matches1, matches2;

  maximal_common_fragments(
      mol1, mol2, frag_noopt, matches1, matches2, 0, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_deg1, matches1, matches2, 0, 0,
      Product::GenerationType::DEG_1,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_uncon, matches1, matches2, 0, 0,
      Product::GenerationType::UNCON,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_uncon_deg1, matches1, matches2, 0, 0,
      Product::GenerationType::UNCON_DEG_1,
      false, false, TIMEOUT_BIG);

  REQUIRE(!frag_noopt.empty());
  REQUIRE(!frag_deg1.empty());
  REQUIRE(!frag_uncon.empty());
  REQUIRE(!frag_uncon_deg1.empty());

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag_noopt.begin(), frag_noopt.end(), comparator);
  std::sort(frag_deg1.begin(), frag_deg1.end(), comparator);
  std::sort(frag_uncon.begin(), frag_uncon.end(), comparator);
  std::sort(frag_uncon_deg1.begin(), frag_uncon_deg1.end(), comparator);

  REQUIRE(frag_noopt[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_deg1[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon[0]->get_atom_count() == mol1.get_atom_count());
  REQUIRE(frag_uncon_deg1[0]->get_atom_count() == mol1.get_atom_count());

}

TEST_CASE("mcf_subisomorphic_graphs_2", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHYL);
  mol2.read_lgf(ETHANE_2);

  FragmentVector frag_noopt, frag_deg1, frag_uncon, frag_uncon_deg1;
  MatchVector matches1, matches2;

  maximal_common_fragments(
      mol1, mol2, frag_noopt, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::DEG_1,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_uncon, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON,
      false, false, TIMEOUT_BIG);
  maximal_common_fragments(
      mol1, mol2, frag_uncon_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON_DEG_1,
      false, false, TIMEOUT_BIG);

  REQUIRE(!frag_noopt.empty());
  REQUIRE(!frag_deg1.empty());
  REQUIRE(!frag_uncon.empty());
  REQUIRE(!frag_uncon_deg1.empty());

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag_noopt.begin(), frag_noopt.end(), comparator);
  std::sort(frag_deg1.begin(), frag_deg1.end(), comparator);
  std::sort(frag_uncon.begin(), frag_uncon.end(), comparator);
  std::sort(frag_uncon_deg1.begin(), frag_uncon_deg1.end(), comparator);

  REQUIRE(frag_noopt[0]->get_atom_count() == 2);
  REQUIRE(frag_deg1[0]->get_atom_count() == 2);
  REQUIRE(frag_uncon[0]->get_atom_count() == 2);
  REQUIRE(frag_uncon_deg1[0]->get_atom_count() == 2);

}

TEST_CASE("mcf_graphs_no_match", "[algo]") {

  Molecule mol1, mol2;
  LGFIOConfig config("label", "atomType");

  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(GENERIC_1, config);

  FragmentVector frag_noopt, frag_deg1, frag_uncon, frag_uncon_deg1;
  MatchVector matches1, matches2;

  maximal_common_fragments(
      mol1, mol2, frag_noopt, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);
  maximal_common_fragments(
      mol1, mol2, frag_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::DEG_1,
      false, false, TIMEOUT);
  maximal_common_fragments(
      mol1, mol2, frag_uncon, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON,
      false, false, TIMEOUT);
  maximal_common_fragments(
      mol1, mol2, frag_uncon_deg1, matches1, matches2, 1, 0,
      Product::GenerationType::UNCON_DEG_1,
      false, false, TIMEOUT);

  REQUIRE(frag_noopt.empty());
  REQUIRE(frag_deg1.empty());
  REQUIRE(frag_uncon.empty());
  REQUIRE(frag_uncon_deg1.empty());

}

TEST_CASE("atomic_fragments", "[algo]") {

  Molecule mol;
  mol.read_lgf(ETHANE_1);

  FragmentVector frag0, frag1, frag2;
  MatchVector matches;

  atomic_fragments(mol, frag0, matches, 0);
  atomic_fragments(mol, frag1, matches, 1);
  atomic_fragments(mol, frag2, matches, 2);

  REQUIRE(frag0.size() == 8);
  REQUIRE(frag1.size() == 8);
  REQUIRE(frag2.size() == 8);

}
