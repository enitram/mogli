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
#include "util/packing.h"

using namespace Catch::Matchers;
using namespace mogli;

TEST_CASE("hash_canonization", "[packing]") {

  Molecule mol1, mol2, mol3;
  LGFIOConfig config("label", "atomType");
  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);
  mol3.read_lgf(GENERIC_1, config);

  Canonization _c1(mol1), _c2(mol2), _c3(mol3);

  auto c1 = hash_canonization(_c1);
  auto c2 = hash_canonization(_c2);
  auto c3 = hash_canonization(_c3);

  REQUIRE(c1 != c2);
  REQUIRE(c1 != c3);

}

TEST_CASE("hash_fcanonization", "[packing]") {

  Molecule mol1, mol2, mol3, mol4;
  LGFIOConfig config("label", "atomType");
  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_1);
  mol3.read_lgf(ETHANE_2);
  mol4.read_lgf(ETHYL, config);

  FragmentVector frag1, frag2, frag3;
  MatchVector matches1, matches2;

  auto t1 = maximal_common_fragments(
      mol1, mol2, frag1, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);
  auto t2 = maximal_common_fragments(
      mol1, mol3, frag2, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);
  auto t3 = maximal_common_fragments(
      mol1, mol4, frag3, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);

  REQUIRE(t1);
  REQUIRE(t2);
  REQUIRE(t3);

  REQUIRE(!frag1.empty());
  REQUIRE(!frag2.empty());
  REQUIRE(!frag3.empty());

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag1.begin(), frag1.end(), comparator);
  std::sort(frag2.begin(), frag2.end(), comparator);
  std::sort(frag3.begin(), frag3.end(), comparator);

  FragmentCanonization _f1(*frag1[0]), _f2(*frag2[0]), _f3(*frag3[0]);

  auto f1 = hash_fcanonization(_f1);
  auto f2 = hash_fcanonization(_f2);
  auto f3 = hash_fcanonization(_f3);

  REQUIRE(f1 == f2);
  REQUIRE(f1 != f3);

}

TEST_CASE("pack_canonization", "[packing]") {

  Molecule mol;
  mol.read_lgf(ETHANE_1);

  Canonization c1(mol), c2;

  unpack_object(pack_object(c1), c2);

  REQUIRE_THAT(c1.get_colors(), Equals(c2.get_colors()));
  REQUIRE_THAT(c1.get_canonization(), Equals(c2.get_canonization()));
  REQUIRE_THAT(c1.get_node_order(), Equals(c2.get_node_order()));

}

TEST_CASE("pack_fcanonization", "[packing]") {

  Molecule mol1, mol2;
  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);

  FragmentVector frag;
  MatchVector matches1, matches2;

  auto t1 = maximal_common_fragments(
      mol1, mol2, frag, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);

  REQUIRE(t1);

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag.begin(), frag.end(), comparator);

  FragmentCanonization f1(*frag[0]), f2;

  unpack_object(pack_object(f1), f2);

  REQUIRE_THAT(f1.get_colors(), Equals(f2.get_colors()));
  REQUIRE_THAT(f1.get_canonization(), Equals(f2.get_canonization()));
  REQUIRE_THAT(f1.get_node_order(), Equals(f2.get_node_order()));
  REQUIRE_THAT(f1.get_core_nodes(), Equals(f2.get_core_nodes()));

}

TEST_CASE("pack_fragment", "[packing]") {

  Molecule mol1, mol2;
  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);

  FragmentVector frag;
  MatchVector matches1, matches2;

  auto t1 = maximal_common_fragments(
      mol1, mol2, frag, matches1, matches2, 1, 0,
      Product::GenerationType::NO_OPT,
      false, false, TIMEOUT);

  REQUIRE(t1);

  auto comparator = [](const std::shared_ptr<Fragment> & a, const std::shared_ptr<Fragment> & b) {
    return a->get_atom_count() > b->get_atom_count();
  };

  std::sort(frag.begin(), frag.end(), comparator);

  Fragment & f1 = *frag[0];
  Fragment f2;
  unpack_object(pack_object(f1), f2);

  REQUIRE(f1.get_atom_count() == f2.get_atom_count());
  REQUIRE(f1.get_core_atom_count() == f2.get_core_atom_count());

  std::vector<std::pair<int, int>> edges1, edges2;
  for (EdgeIt e = f1.get_edge_iter(); e != lemon::INVALID; ++e) {
    int id1 = f1.get_id(f1.get_u(e));
    int id2 = f1.get_id(f1.get_v(e));
    if (id1 < id2) {
      edges1.emplace_back(id1, id2);
    } else {
      edges1.emplace_back(id2, id1);
    }
  }
  for (EdgeIt e = mol2.get_edge_iter(); e != lemon::INVALID; ++e) {
    int id1 = f2.get_id(f2.get_u(e));
    int id2 = f2.get_id(f2.get_v(e));
    if (id1 < id2) {
      edges2.emplace_back(id1, id2);
    } else {
      edges2.emplace_back(id2, id1);
    }
  }

  std::sort(edges1.begin(), edges1.end());
  std::sort(edges2.begin(), edges2.end());

  REQUIRE_THAT(edges1, Equals(edges2));

  for (NodeIt v1 = f1.get_node_iter(); v1 != lemon::INVALID; ++v1) {
    auto v2 = f2.get_node_by_id(f1.get_id(v1));
    REQUIRE(f1.is_core(v1) == f2.is_core(v2));
  }

}

TEST_CASE("pack_molecule", "[packing]") {
  
  Molecule mol1, mol2;
  LGFIOConfig config("label", "atomType");
  config.add_string_node_prop("label2").add_bool_node_prop("isC").add_int_node_prop("iNum").add_double_node_prop("dNum");
  
  mol1.read_lgf(ETHANE_PROPS, config);
  
  unpack_object(pack_object(mol1), mol2);

  REQUIRE(mol1.get_atom_count() == mol2.get_atom_count());

  std::vector<std::pair<int, int>> edges1, edges2;
  for (EdgeIt e = mol1.get_edge_iter(); e != lemon::INVALID; ++e) {
    int id1 = mol1.get_id(mol1.get_u(e));
    int id2 = mol1.get_id(mol1.get_v(e));
    if (id1 < id2) {
      edges1.emplace_back(id1, id2);
    } else {
      edges1.emplace_back(id2, id1);
    }
  }
  for (EdgeIt e = mol2.get_edge_iter(); e != lemon::INVALID; ++e) {
    int id1 = mol2.get_id(mol2.get_u(e));
    int id2 = mol2.get_id(mol2.get_v(e));
    if (id1 < id2) {
      edges2.emplace_back(id1, id2);
    } else {
      edges2.emplace_back(id2, id1);
    }
  }

  std::sort(edges1.begin(), edges1.end());
  std::sort(edges2.begin(), edges2.end());

  REQUIRE_THAT(edges1, Equals(edges2));

  StringVector props;
  mol2.get_properties(props);
  REQUIRE_THAT(props, VectorContains(std::string("label2")));
  REQUIRE_THAT(props, VectorContains(std::string("isC")));
  REQUIRE_THAT(props, VectorContains(std::string("iNum")));
  REQUIRE_THAT(props, VectorContains(std::string("dNum")));

  for (NodeIt v1 = mol1.get_node_iter(); v1 != lemon::INVALID; ++v1) {
    auto v2 = mol2.get_node_by_id(mol1.get_id(v1));
    REQUIRE(mol1.get_property(v1, "label2") == mol2.get_property(v2, "label2"));
    REQUIRE(mol1.get_property(v1, "isC") == mol2.get_property(v2, "isC"));
    REQUIRE(mol1.get_property(v1, "iNum") == mol2.get_property(v2, "iNum"));
    REQUIRE(mol1.get_property(v1, "dNum") == mol2.get_property(v2, "dNum"));
  }

}

TEST_CASE("pack_match", "[packing]") {

  Match m1, m2;

  m1.add_frag_to_mol(0, 1);
  m1.add_frag_to_mol(1, 2);

  unpack_object(pack_object(m1), m2);

  IntVector i1, i2;
  m1.get_atom_ids(i1);
  m2.get_atom_ids(i2);

  REQUIRE_THAT(i1, Equals(i2));
  REQUIRE(m1.frag_to_mol(0) == m2.frag_to_mol(0));
  REQUIRE(m1.frag_to_mol(1) == m2.frag_to_mol(1));
  REQUIRE(m1.frag_to_mol(2) == m2.frag_to_mol(2));

}
