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

#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "mcf.h"
#include "util/packing.h"


namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(mogli::BoolVector);
PYBIND11_MAKE_OPAQUE(mogli::ShortVector);
PYBIND11_MAKE_OPAQUE(mogli::IntVector);
PYBIND11_MAKE_OPAQUE(mogli::LongVector);
PYBIND11_MAKE_OPAQUE(mogli::StringVector);
PYBIND11_MAKE_OPAQUE(mogli::FragmentVector);
PYBIND11_MAKE_OPAQUE(mogli::MatchVector);
PYBIND11_MAKE_OPAQUE(mogli::SharedPtrVector<mogli::Molecule>::type)
PYBIND11_MAKE_OPAQUE(mogli::NodeVector);

inline py::object pass_through(py::object const& o) { return o; }

template <typename T>
static std::shared_ptr<T> deserialize(std::string str) {
  auto obj = std::make_shared<T>();
  mogli::unpack_object(str, *obj);
  return obj;
}

template <typename T, typename E>
struct ItWrapper {

  static E next(T& o) {
    E next = o;
    if (next == lemon::INVALID) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      throw py::error_already_set();
    }
    ++o;
    return next;
  }

  static void wrap(py::handle handle, const char * name) {
    py::class_<T>(handle, name)
        .def("next", next)
        .def("__next__", next)
        .def("__iter__", pass_through);
  }
};

PYBIND11_MODULE(mogli, m) {

  // gaph iterators

  ItWrapper<mogli::EdgeIt, const mogli::Edge>::wrap(m, "EdgeIt");
  ItWrapper<mogli::IncEdgeIt, const mogli::Edge>::wrap(m, "IncEdgeIt");
  ItWrapper<mogli::NodeIt, const mogli::Node>::wrap(m, "NodeIt");

  // STL containers

  py::bind_vector<mogli::BoolVector>(m, "BoolVector");
  py::bind_vector<mogli::ShortVector>(m, "ShortVector", py::buffer_protocol());
  py::bind_vector<mogli::IntVector>(m, "IntVector", py::buffer_protocol());
  py::bind_vector<mogli::LongVector>(m, "LongVector", py::buffer_protocol());
  py::bind_vector<mogli::StringVector>(m, "StringVector");
  py::bind_vector<mogli::FragmentVector>(m, "FragmentVector");
  py::bind_vector<mogli::MatchVector>(m, "MatchVector");
  py::bind_vector<mogli::SharedPtrVector<mogli::Molecule>::type>(m, "MoleculeVector");
  py::bind_vector<mogli::NodeVector>(m, "NodeVector");

  // enums

  py::enum_<mogli::Product::GenerationType>(m, "GenerationType")
      .value("NO_OPT", mogli::Product::GenerationType::NO_OPT, "No optimization")
      .value("UNCON", mogli::Product::GenerationType::UNCON, "Graph partitioning rule")
      .value("DEG_1", mogli::Product::GenerationType::DEG_1, "Degree-1 neighbors rule")
      .value("UNCON_DEG_1", mogli::Product::GenerationType::UNCON_DEG_1, "Degree-1 neighbors and graph partitioning rules");

  // classes

  py::class_<mogli::Canonization>(m, "Canonization")
      .def(py::init<const mogli::Molecule&>())
      .def("get_colors", &mogli::Canonization::get_colors, py::return_value_policy::reference_internal, "Get element numbers of the ordered atoms")
      .def("get_canonization", &mogli::Canonization::get_canonization, py::return_value_policy::reference_internal, "Get canonization")
      .def("get_node_order", &mogli::Canonization::get_node_order, py::return_value_policy::reference_internal, "Get ordered atom ids")
      .def("is_isomorphic", &mogli::Canonization::is_isomorphic, "Isomorphism test");

  py::class_<mogli::Edge>(m, "Edge")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self);

  py::class_<mogli::FragmentCanonization, mogli::Canonization>(m, "FragmentCanonization")
      .def(py::init<const mogli::Fragment&>())
      .def("get_core_nodes", &mogli::FragmentCanonization::get_core_nodes, py::return_value_policy::reference_internal, "Get atom ids of the ordered core atoms")
      .def("is_isomorphic", &mogli::FragmentCanonization::is_isomorphic, "Isomorphism test");

  py::class_<mogli::LGFIOConfig>(m, "LGFIOConfig")
      .def(py::init<std::string, std::string>())
      .def("add_bool_node_prop", &mogli::LGFIOConfig::add_bool_node_prop, py::return_value_policy::reference_internal, "Add a node property that will be cast to a bool")
      .def("add_int_node_prop", &mogli::LGFIOConfig::add_int_node_prop, py::return_value_policy::reference_internal, "Add a node property that will be cast to an int")
      .def("add_double_node_prop", &mogli::LGFIOConfig::add_double_node_prop, py::return_value_policy::reference_internal, "Add a node property that will be cast to a double")
      .def("add_string_node_prop", &mogli::LGFIOConfig::add_string_node_prop, py::return_value_policy::reference_internal, "Add a node property that will be cast to a str");

  py::class_<mogli::Match, std::shared_ptr<mogli::Match>>(m, "Match")
      .def("frag_to_mol", &mogli::Match::frag_to_mol)
      .def("merged_frag_to_mol",
          [](const mogli::Match &self, const int id) {
            mogli::IntVector ids;
            self.merged_frag_to_mol(id, ids);
            return ids;
          },
          py::return_value_policy::move
          )
      .def("get_atom_ids",
          [](const mogli::Match & self) {
            mogli::IntVector ids;
            self.get_atom_ids(ids);
            return ids;
          },
          py::return_value_policy::move
          )
      .def("map_ids", &mogli::Match::map_ids);

  py::class_<mogli::Molecule, std::shared_ptr<mogli::Molecule>> molecule(m, "Molecule");
  molecule.def(py::init<>())
      .def(py::init<mogli::PeriodicTable&>())
      .def("add_atom", py::overload_cast<std::string>(&mogli::Molecule::add_atom), "Add atom with element string")
      .def("add_atom", py::overload_cast<int, std::string>(&mogli::Molecule::add_atom), "Add atom with element string and node id")
      .def("add_atom", py::overload_cast<unsigned short>(&mogli::Molecule::add_atom), "Add atom with element number")
      .def("add_atom", py::overload_cast<int, unsigned short>(&mogli::Molecule::add_atom), "Add atom with element number and node id")
      .def("add_edge", &mogli::Molecule::add_edge, "Add bond")
      .def("get_atom_count", &mogli::Molecule::get_atom_count, "Get number of atoms")
      .def("get_color", &mogli::Molecule::get_color, "Get element number of this atom")
      .def("get_connected_components",
          [](mogli::Molecule &self) {
            mogli::SharedPtrVector<mogli::Molecule>::type components;
            self.get_connected_components(components);
            return components;
          },
          py::return_value_policy::move,
          "Get connected components of the molecular graph")
      .def("get_edge_iter", &mogli::Molecule::get_edge_iter, "Get bond iterator")
      .def("get_element", &mogli::Molecule::get_element, "Get element of this atom")
      .def("get_id", &mogli::Molecule::get_id, "Get id of this atom")
      .def("get_inc_edge_iter", &mogli::Molecule::get_inc_edge_iter, "Get incident bonds iterator")
      .def("get_node_by_id", &mogli::Molecule::get_node_by_id, "Get atom by id")
      .def("get_node_iter", &mogli::Molecule::get_node_iter, "Get atom iterator")
      .def("get_opposite_node", &mogli::Molecule::get_opposite_node, "Get atom on the opposite side of the bond")
      .def("get_properties",
          [](const mogli::Molecule &self) {
            mogli::StringVector properties;
            self.get_properties(properties);
            return properties;
          },
          py::return_value_policy::move,
          "Get all property names")
      .def("get_property", &mogli::Molecule::get_property, "Get atom property")
      .def("get_u", &mogli::Molecule::get_u, "Get first atom of the bond")
      .def("get_v", &mogli::Molecule::get_v, "Get second atom of the bond")
      .def("has_node_with_id", &mogli::Molecule::has_node_with_id, "Tests if atom with this id exists")
      .def("is_connected", &mogli::Molecule::is_connected, "Test if molecular graph is connected")
      .def("is_isomorphic", &mogli::Molecule::is_isomorphic, "Isomorphisms test")
      .def("print_dot", py::overload_cast<>(&mogli::Molecule::print_dot, py::const_), "Print dot (graphviz file)")
      .def("print_dot", py::overload_cast<const mogli::StringVector&>(&mogli::Molecule::print_dot, py::const_), "Print dot (graphviz) file with atom properties")
      .def("read_lgf", py::overload_cast<const std::string &>(&mogli::Molecule::read_lgf), "Read default formatted LGF file")
      .def("read_lgf", py::overload_cast<const std::string &, const mogli::LGFIOConfig&>(&mogli::Molecule::read_lgf), "Read LGF file")
      .def("set_property", &mogli::Molecule::set_property, "Set atom property")
      .def("split", &mogli::Molecule::split, "Balanced split of the molecule")
      .def("write_lgf", py::overload_cast<>(&mogli::Molecule::write_lgf), "Write default formatted LGF file")
      .def("write_lgf", py::overload_cast<const mogli::LGFIOConfig&>(&mogli::Molecule::write_lgf), "Write LGF file");

  py::class_<mogli::Fragment, std::shared_ptr<mogli::Fragment>>(m, "Fragment", molecule)
      .def("get_core_atom_count", &mogli::Fragment::get_core_atom_count, "Get number of core atoms")
      .def("is_core", &mogli::Fragment::is_core, "Is this node a core atom?")
      .def("print_dot", py::overload_cast<>(&mogli::Fragment::print_dot,  py::const_), "Print dot (graphviz) file");

  py::class_<mogli::PeriodicTable>(m, "PeriodicTable")
      .def(py::init<>())
      .def("add", &mogli::PeriodicTable::add, py::return_value_policy::reference_internal, "Add a new element")
      .def("add_uncolored", &mogli::PeriodicTable::add_uncolored, py::return_value_policy::reference_internal, "Add a new element with an associated color string (for dot printing)");

  py::class_<mogli::Node>(m, "Node")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self);

  // methods

  m.def("atomic_fragments",
      [](mogli::Molecule &mol, int shell) {
        mogli::FragmentVector fragments;
        mogli::MatchVector matches;
        mogli::atomic_fragments(mol, fragments, matches, shell);
        return py::make_tuple(fragments, matches);
      },
      py::return_value_policy::move,
      "Get all single-atom fragments");

  m.def("maximal_common_fragments",
      [](mogli::Molecule &mol1, mogli::Molecule &mol2,
          int shell, unsigned int min_core_size, unsigned int max_core_size,
          mogli::Product::GenerationType prod_gen, bool maximum, int timeout_seconds){
        mogli::FragmentVector fragments;
        mogli::MatchVector matches_mol1, matches_mol2;
        mogli::maximal_common_fragments(mol1, mol2, fragments, matches_mol1, matches_mol2,
            shell, min_core_size, max_core_size, prod_gen, false, maximum, timeout_seconds);
        return py::make_tuple(fragments, matches_mol1, matches_mol2);
      },
      py::return_value_policy::move,
      "Compute maximal common fragments");

  m.def("maximal_common_fragments",
      [](mogli::Molecule &mol1, mogli::Molecule &mol2,
         int shell, unsigned int min_core_size,
         mogli::Product::GenerationType prod_gen, bool maximum, int timeout_seconds){
        mogli::FragmentVector fragments;
        mogli::MatchVector matches_mol1, matches_mol2;
        mogli::maximal_common_fragments(mol1, mol2, fragments, matches_mol1, matches_mol2,
                                        shell, min_core_size, prod_gen, false, maximum, timeout_seconds);
        return py::make_tuple(fragments, matches_mol1, matches_mol2);
      },
      py::return_value_policy::move,
      "Compute maximal common fragments");

  // hashing methods

  m.def("hash_canonization", &mogli::hash_canonization, "Get molecule canonization hash");

  m.def("hash_fcanonization", &mogli::hash_fcanonization, "Get fragment canonization hash");

  // serializing methods

  m.def("pack_canonization", py::overload_cast<const mogli::Canonization&>(&mogli::pack_object<mogli::Canonization>), "Serialize molecule canonization");

  m.def("pack_canonization", py::overload_cast<const std::shared_ptr<mogli::Canonization>&>(&mogli::pack_object<mogli::Canonization>), "Serialize molecule canonization");

  m.def("pack_fcanonization", py::overload_cast<const mogli::FragmentCanonization&>(&mogli::pack_object<mogli::FragmentCanonization>), "Serialize fragment canonization");

  m.def("pack_fcanonization", py::overload_cast<const std::shared_ptr<mogli::FragmentCanonization>&>(&mogli::pack_object<mogli::FragmentCanonization>), "Serialize fragment canonization");

  m.def("pack_fragment", py::overload_cast<const mogli::Fragment&>(&mogli::pack_object<mogli::Fragment>), "Serialize fragment");

  m.def("pack_fragment", py::overload_cast<const std::shared_ptr<mogli::Fragment>&>(&mogli::pack_object<mogli::Fragment>), "Serialize fragment");

  m.def("pack_match", py::overload_cast<const mogli::Match&>(&mogli::pack_object<mogli::Match>), "Serialize match object");

  m.def("pack_match", py::overload_cast<const std::shared_ptr<mogli::Match>&>(&mogli::pack_object<mogli::Match>), "Serialize match object");

  m.def("pack_molecule", py::overload_cast<const mogli::Molecule&>(&mogli::pack_object<mogli::Molecule>), "Serialize molecule");

  m.def("pack_molecule", py::overload_cast<const std::shared_ptr<mogli::Molecule>&>(&mogli::pack_object<mogli::Molecule>), "Serialize molecule");

  // de-serializing methods

  m.def("unpack_canonization", &deserialize<mogli::Canonization>, "Deserialize molecule canonization");

  m.def("unpack_fcanonization", &deserialize<mogli::FragmentCanonization>, "Deserizalize fragment canonization");

  m.def("unpack_fragment", &deserialize<mogli::Fragment>, "Deserizalize fragment");

  m.def("unpack_match", &deserialize<mogli::Match>, "Deserialize match object");

  m.def("unpack_molecule", &deserialize<mogli::Molecule>, "Deserialize molecule");

}