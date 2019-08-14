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
#include "fcanonization.h"
#include "mcf.h"
#include "test_fixtures.h"

using namespace Catch::Matchers;
using namespace mogli;

TEST_CASE("canonization", "[algo]") {

  Molecule mol1, mol2, mol3;
  LGFIOConfig config("label", "atomType");
  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);
  mol3.read_lgf(GENERIC_1, config);

  Canonization c1(mol1), c2(mol2), c3(mol3);

  REQUIRE_THAT(c1.get_canonization(), Equals(c2.get_canonization()));
  REQUIRE_THAT(c1.get_colors(), Equals(c2.get_colors()));
  REQUIRE_THAT(c1.get_canonization(), !Equals(c3.get_canonization()));
  REQUIRE_THAT(c1.get_colors(), !Equals(c3.get_colors()));

}

TEST_CASE("fcanonization", "[algo]") {

  Molecule mol1, mol2, mol3, mol4;
  LGFIOConfig config("label", "atomType");
  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_1);
  mol3.read_lgf(ETHANE_1);
  mol4.read_lgf(ETHYL, config);

  FragmentVector frag1, frag2, frag3;
  MatchVector matches1, matches2;

  auto t1 = maximal_common_fragments(
      mol1, mol2, frag1, matches1, matches2, 1, TIMEOUT,Product::GenerationType::NO_OPT);
  auto t2 = maximal_common_fragments(
      mol1, mol3, frag2, matches1, matches2, 1, TIMEOUT,Product::GenerationType::NO_OPT);
  auto t3 = maximal_common_fragments(
      mol1, mol4, frag3, matches1, matches2, 1, TIMEOUT,Product::GenerationType::NO_OPT);

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

  FragmentCanonization f1(*frag1[0]), f2(*frag2[0]), f3(*frag3[0]);

  REQUIRE_THAT(f1.get_canonization(), Equals(f2.get_canonization()));
  REQUIRE_THAT(f1.get_colors(), Equals(f2.get_colors()));
  REQUIRE_THAT(f1.get_canonization(), !Equals(f3.get_canonization()));
  REQUIRE_THAT(f1.get_colors(), !Equals(f3.get_colors()));

}

TEST_CASE("isomorphism", "[algo]") {

  Molecule mol1, mol2, mol3;

  mol1.read_lgf(ETHANE_1);
  mol2.read_lgf(ETHANE_2);
  mol3.read_lgf(ETHYL);

  REQUIRE(mol1.is_isomorphic(mol2));
  REQUIRE(!mol1.is_isomorphic(mol3));

}

TEST_CASE("subgraph_isomorphism", "[algo]") {

  Molecule mol1, mol2;

  mol1.read_lgf(ETHYL);
  mol2.read_lgf(ETHANE_2);

  IntToIntMap iso;

  REQUIRE(are_subgraph_isomorphic(mol1, mol2, iso));
  REQUIRE(!are_subgraph_isomorphic(mol2, mol1, iso));

}
