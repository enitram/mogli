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

using namespace Catch::Matchers;
using namespace mogli;

TEST_CASE("read_lgf", "[io]") {

  Molecule mol;
  mol.read_lgf(ETHANE_1);

  REQUIRE(mol.get_atom_count() == 8);

  auto c1 = mol.get_node_by_id(0);
  auto c2 = mol.get_node_by_id(1);

  StringVector props;
  mol.get_properties(props);

  REQUIRE_THAT(props, VectorContains(std::string("label2")));
  REQUIRE(std::get<std::string>(mol.get_property(c1, "label2")) == "C1");
  REQUIRE(std::get<std::string>(mol.get_property(c2, "label2")) == "C2");

  std::vector<std::pair<int, int>> c1_edges;
  c1_edges.emplace_back(0, 1);
  c1_edges.emplace_back(0, 2);
  c1_edges.emplace_back(0, 3);
  c1_edges.emplace_back(0, 4);

  std::vector<std::pair<int, int>> c2_edges;
  c2_edges.emplace_back(0, 1);
  c2_edges.emplace_back(1, 5);
  c2_edges.emplace_back(1, 6);
  c2_edges.emplace_back(1, 7);

  std::vector<std::pair<int, int>> edges(c1_edges);
  std::copy(c2_edges.begin(), c2_edges.end(), std::inserter(edges, edges.end()));

  for (auto e = mol.get_inc_edge_iter(c1); e != lemon::INVALID; ++e) {
    auto id1 = mol.get_id(mol.get_u(e));
    auto id2 = mol.get_id(mol.get_v(e));
    auto edge = (id1 < id2) ? std::make_pair(id1, id2) : std::make_pair(id2, id1);
    REQUIRE_THAT(c1_edges, VectorContains(edge));
  }

  for (auto e = mol.get_inc_edge_iter(c2); e != lemon::INVALID; ++e) {
    auto id1 = mol.get_id(mol.get_u(e));
    auto id2 = mol.get_id(mol.get_v(e));
    auto edge = (id1 < id2) ? std::make_pair(id1, id2) : std::make_pair(id2, id1);
    REQUIRE_THAT(c2_edges, VectorContains(edge));
  }

  for (auto e = mol.get_edge_iter(); e != lemon::INVALID; ++e) {
    auto id1 = mol.get_id(mol.get_u(e));
    auto id2 = mol.get_id(mol.get_v(e));
    auto edge = (id1 < id2) ? std::make_pair(id1, id2) : std::make_pair(id2, id1);
    REQUIRE_THAT(edges, VectorContains(edge));
  }

}

TEST_CASE("read_lgf_properties", "[io]") {

  LGFIOConfig config("label", "atomType");
  config.add_string_node_prop("label2").add_bool_node_prop("isC").add_int_node_prop("iNum").add_double_node_prop("dNum");

  Molecule mol;
  mol.read_lgf(ETHANE_PROPS, config);

  REQUIRE(mol.get_atom_count() == 8);

  StringVector props;
  mol.get_properties(props);

  REQUIRE(props.size() == 4);

  REQUIRE_THAT(props, VectorContains(std::string("label2")));
  REQUIRE_THAT(props, VectorContains(std::string("isC")));
  REQUIRE_THAT(props, VectorContains(std::string("iNum")));
  REQUIRE_THAT(props, VectorContains(std::string("dNum")));

  auto v = mol.get_node_by_id(0);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "C1");
  REQUIRE(std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 1);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 1.0);

  v = mol.get_node_by_id(1);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "C2");
  REQUIRE(std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 2);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 2.0);

  v = mol.get_node_by_id(2);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "H3");
  REQUIRE(!std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 3);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 3.0);

  v = mol.get_node_by_id(3);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "H4");
  REQUIRE(!std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 4);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 4.0);

  v = mol.get_node_by_id(4);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "H5");
  REQUIRE(!std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 5);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 5.0);

  v = mol.get_node_by_id(5);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "H6");
  REQUIRE(!std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 6);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 6.0);

  v = mol.get_node_by_id(6);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "H7");
  REQUIRE(!std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 7);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 7.0);

  v = mol.get_node_by_id(7);
  REQUIRE(std::get<std::string>(mol.get_property(v, "label2")) == "H8");
  REQUIRE(!std::get<bool>(mol.get_property(v, "isC")));
  REQUIRE(std::get<int>(mol.get_property(v, "iNum")) == 8);
  REQUIRE(std::get<double>(mol.get_property(v, "dNum")) == 8.0);

}
