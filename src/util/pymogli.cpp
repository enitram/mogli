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
using namespace py::literals;

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
static py::bytes serialize(const T & obj) {
  return py::bytes(mogli::pack_object(obj));
}

template <typename T>
static std::shared_ptr<T> deserialize_ptr(const std::string & str) {
  auto obj = std::make_shared<T>();
  mogli::unpack_object(str, *obj);
  return obj;
}

template <typename T>
static T deserialize(const std::string & str) {
  T obj;
  mogli::unpack_object(str, obj);
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
  py::bind_map<mogli::IntToIntMap>(m, "IntToIntMap");

  // enums

  py::enum_<mogli::Product::GenerationType>(m, "GenerationType")
      .value("NO_OPT", mogli::Product::GenerationType::NO_OPT, "No optimization")
      .value("UNCON", mogli::Product::GenerationType::UNCON, "Graph partitioning rule")
      .value("DEG_1", mogli::Product::GenerationType::DEG_1, "Degree-1 neighbors rule")
      .value("UNCON_DEG_1", mogli::Product::GenerationType::UNCON_DEG_1, "Degree-1 neighbors and graph partitioning rules");

  // classes

  py::class_<mogli::Canonization> canon(m, "Canonization");

  canon.doc() = "Canonical representation of a molecular graph.";

  canon.def(
          py::init<const mogli::Molecule&>(),
          "mol"_a,
          R"(
          Create a canonical representation of a molecular graph.

          Args:
              mol (Molecule): Molecular graph.
          )")
      .def(
          "get_colors",
          &mogli::Canonization::get_colors,
          py::return_value_policy::reference_internal,
          R"(
          Returns the element numbers of the atoms in canonical order.

          Returns:
              ShortVector. Element numbers.
          )")
      .def("get_canonization",
          &mogli::Canonization::get_canonization,
          py::return_value_policy::reference_internal,
          R"(
          Returns the canonical representations of the atoms.

          Returns:
              LongVector. Canonical representations.
          )")
      .def("get_node_order",
          &mogli::Canonization::get_node_order,
          py::return_value_policy::reference_internal,
          R"(
          Returns the atom IDs in canonical order.

          Returns:
              ShortVector. Atom IDs.
          )")
      .def("is_isomorphic",
          &mogli::Canonization::is_isomorphic,
          "other"_a,
          R"(
          Isomorphism test.

          Args:
              other (Canonization): Other canonization.

          Returns:
              bool. True, if isomorphic to other canonization, false otherwise.
          )");

  py::class_<mogli::Edge>(m, "Edge")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self);

  py::class_<mogli::FragmentCanonization, mogli::Canonization> fcanon(m, "FragmentCanonization");

  fcanon.doc() = "Canonical representation of a molecular fragment.";

  fcanon.def(
          py::init<const mogli::Fragment&>(),
          "fragment"_a,
          R"(
          Create a canonical representation of a molecular fragment.

          Args:
              fragment (Fragment): Molecular fragment.
          )")
      .def(
          "get_core_nodes",
          &mogli::FragmentCanonization::get_core_nodes,
          py::return_value_policy::reference_internal,
          R"(
          Returns a bool vector indicating the core atoms in canonical order.

          Returns:
              BoolVector. Core atoms.
          )")
      .def(
          "is_isomorphic",
          &mogli::FragmentCanonization::is_isomorphic,
          "other"_a,
          R"(
          Isomorphism test.

          Args:
              other (Canonization): Other canonization.

          Returns:
              bool. True, if isomorphic to other canonization, false otherwise.

          )");

  py::class_<mogli::LGFIOConfig> conf(m, "LGFIOConfig");

  conf.doc() = "LGF-formatter for reading and writing molecular graphs.";

  conf.def(
          py::init<std::string, std::string>(),
          "id_property"_a, "color_property"_a,
          R"(
          Initialize an LGF formatter.

          Args:
              id_property (str):    Name of the atom ID column.
              color_property (str): Name of the element number column.
          )")
      .def(
          "add_bool_node_prop",
          &mogli::LGFIOConfig::add_bool_node_prop,
          "property"_a,
          py::return_value_policy::reference_internal,
          R"(
          Add a bool atom property column.

          Args:
              property (str): Column name.

          Returns:
             LGFIOConfig. Updated LGF formatter.
          )")
      .def(
          "add_int_node_prop",
          &mogli::LGFIOConfig::add_int_node_prop,
          "property"_a,
          py::return_value_policy::reference_internal,
          R"(
          Add an int atom property column.

          Args:
              property (str): Column name.

          Returns:
             LGFIOConfig. Updated LGF formatter.
          )")
      .def(
          "add_double_node_prop",
          &mogli::LGFIOConfig::add_double_node_prop,
          "property"_a,
          py::return_value_policy::reference_internal,
          R"(
          Add a double atom property column.

          Args:
              property (str): Column name.

          Returns:
             LGFIOConfig. Updated LGF formatter.
          )")
      .def("add_string_node_prop",
          &mogli::LGFIOConfig::add_string_node_prop,
          "property"_a,
          py::return_value_policy::reference_internal,
           R"(
          Add a string atom property column.

          Args:
              property (str): Column name.

          Returns:
             LGFIOConfig. Updated LGF formatter.
          )");

  py::class_<mogli::Match, std::shared_ptr<mogli::Match>> match(m, "Match");

  match.doc() = "Maps fragment atom IDs to molecule atom IDs.";

  match.def(
          py::init<>(),
          R"(
          Empty constructor.
          )")
      .def(
          "add_frag_to_mol",
          &mogli::Match::add_frag_to_mol,
          "from"_a, "to"_a,
          R"(
          Add a new fragment to molecule atom mapping.

          Args:
              from (int): Fragment atom ID.
              to (int):   Molecule atom ID.)")
      .def("add_merged_frag_to_mol",
          &mogli::Match::add_merged_frag_to_mol,
          "ftm"_a,
          R"(
          Add new fragment to molecule atom mappings.

          Args:
              ftm (IntToIntMap): Fragment to molecule atom mappings.
          )")
      .def("frag_to_mol",
          &mogli::Match::frag_to_mol,
          "id"_a,
          R"(
          Map fragment atom ID to molecule atom ID.

          Args:
              id (int): Fragment atom ID.

          Returns:
              int. Molecule atom ID.)")
      .def("merged_frag_to_mol",
          [](const mogli::Match &self, const int id) {
            mogli::IntVector ids;
            self.merged_frag_to_mol(id, ids);
            return ids;
          },
          "id"_a,
          py::return_value_policy::move,
          R"(
          Map fragment atom ID to molecule atom IDs.

          Args:
              id (int): Fragment atom ID.

          Returns:
              IntVector. Molecule atom IDs.
          )")
      .def("get_atom_ids",
          [](const mogli::Match & self) {
            mogli::IntVector ids;
            self.get_atom_ids(ids);
            return ids;
          },
          py::return_value_policy::move,
          R"(
          Returns all mapped molecule IDs.

          Returns:
              IntVector. Mapped molecule IDs.
          )")
      .def("merge",
          &mogli::Match::merge,
          "other"_a, "isomorphism_map"_a,
          R"(
          Merge with another match object.

          Args:
              other (Match):                 Other match object.
              isomorphism_map (IntToIntMap): Isomorphism map.)")
      .def("map_ids",
          &mogli::Match::map_ids,
          "this_canon"_a, "other_canon"_a,
          R"(
          Transform fragment to molecule atom mapping to match another molecular graph.

          Args:
              this_canon (Canonization):  Canonization of the molecular graph this match object maps to.
              other_canon (Canonization): Canonization of another molecular graph.
          )");

  py::class_<mogli::Molecule, std::shared_ptr<mogli::Molecule>> molecule(m, "Molecule");

  molecule.doc() = "Molecular graph.";

  molecule.def(
          py::init<>(),
          R"(
          Empty constructor.
          )")
      .def(
          py::init<mogli::PeriodicTable&>(),
          "periodic_table"_a,
          R"(
          Constructor with custom periodic table.

          Args:
              periodic_table (PeriodicTable): Periodic table.
          )")
      .def("add_atom",
          py::overload_cast<std::string>(&mogli::Molecule::add_atom),
          "element"_a,
          R"(
          Add new atom.

          Args:
              element (str): Element type.

          Returns:
              Node. Atom.
          )")
      .def("add_atom",
          py::overload_cast<int, std::string>(&mogli::Molecule::add_atom),
          "id"_a, "element"_a,
               R"(
          Add new atom.

          Args:
              id (int):      Atom ID.
              element (str): Element type.

          Returns:
              Node. Atom.
          )")
      .def("add_atom",
          py::overload_cast<unsigned short>(&mogli::Molecule::add_atom),
          "color"_a,
           R"(
          Add new atom.

          Args:
              color (int): Element number.

          Returns:
              Node. Atom.
          )")
      .def("add_atom",
          py::overload_cast<int, unsigned short>(&mogli::Molecule::add_atom),
          "id"_a, "color"_a,
           R"(
          Add new atom.

          Args:
              id (int):    Atom ID.
              color (int): Element number.

          Returns:
              Node. Atom.
          )")
      .def("add_edge",
          &mogli::Molecule::add_edge,
          "u"_a, "v"_a,
          R"(
          Add bond.

          Args:
              u (Node): Atom.
              v (Node): Atom.
          )")
      .def("get_atom_count",
          &mogli::Molecule::get_atom_count,
          R"(
          Returns the number of atoms.

          Returns:
              int. Number of atoms.
          )")
      .def("get_color",
          &mogli::Molecule::get_color,
          "node"_a,
          R"(
          Returns the element number of this atom.

          Args:
              node (Node): Atom.

          Returns:
              int. Element number.
          )")
      .def("get_connected_components",
          [](mogli::Molecule &self) {
            mogli::SharedPtrVector<mogli::Molecule>::type components;
            self.get_connected_components(components);
            return components;
          },
          py::return_value_policy::move,
          R"(
          Returns the connected components of the molecular graph.

          Returns:
               components (MoleculeVector): Connected components.
          )")
      .def("get_edge_iter",
          &mogli::Molecule::get_edge_iter,
          R"(
          Returns an iterator over all bonds.

          Returns:
              Iterable[Edge]. Bond iterator.
          )")
      .def("get_element",
          &mogli::Molecule::get_element,
          "node"_a,
          R"(
          Returns the element type of this atom.

          Args:
              node (Node): Atom.

          Returns:
              str. Element type.
          )")
      .def("get_id",
          &mogli::Molecule::get_id,
          "node"_a
          R"(
          Returns the ID of this atom.

          Args:
              node (Node): Atom.

          Returns:
              int. Atom ID.
          )")
      .def("get_inc_edge_iter",
          &mogli::Molecule::get_inc_edge_iter,
          "node"_a,
          R"(
          Returns an iterator over all incident bonds.

          Args:
              node (Node): Atom.

          Returns:
              Iterable[Edge]. Bond iterator.
          )")
      .def("get_node_by_id",
          &mogli::Molecule::get_node_by_id,
          "id"_a,
          R"(
          Returns atom with this ID.

          Args:
              id (int): Atom ID.

          Returns:
              Node. Atom.
          )")
      .def("get_node_iter",
          &mogli::Molecule::get_node_iter,
          R"(
          Returns an iterator over all atoms.

          Returns:
              Iterable[Node]. Atom iterator.
          )")
      .def("get_opposite_node",
          &mogli::Molecule::get_opposite_node,
          "node"_a, "edge"_a,
          R"(
          Returns the atom on the opposite side of the bond.

          Args:
              node (Node): Atom.
              edge (Edge): Bond.

          Returns:
              Node. Opposite atom.
          )")
      .def("get_properties",
          [](const mogli::Molecule &self) {
            mogli::StringVector properties;
            self.get_properties(properties);
            return properties;
          },
          py::return_value_policy::move,
          R"(
          Returns the names of all atom properties.

          Returns:
              StringVector. Property names.
          )")
      .def("get_property",
          &mogli::Molecule::get_property,
          "node"_a, "property"_a,
          R"(
          Returns atom property with this name.

          Args:
              node (Node):    Atom.
              property (str): Property name.

          Returns:
              Union[bool, int, float, str]. Property value.
          )")
      .def("get_u",
          &mogli::Molecule::get_u,
          "edge"_a,
          R"(
          Returns the first atom of this bond.

          Args:
              edge (Edge): Bond.

          Returns:
              Node. First atom.
          )")
      .def("get_v",
          &mogli::Molecule::get_v,
          "edge"_a,
           R"(
          Returns the second atom of this bond.

          Args:
              edge (Edge): Bond.

          Returns:
              Node. Second atom.
          )")
      .def("has_node_with_id",
          &mogli::Molecule::has_node_with_id,
          "id"_a,
          R"(
          Test if atom with this ID exists.

          Args:
              id (int): Atom ID.

          Returns:
              bool. True, if atom with this ID exists, false otherwise.
          )")
      .def("is_connected",
          &mogli::Molecule::is_connected,
          R"(
          Test if molecular graph is connected.

          Returns:
              bool. True, if molecular graph is connected, false otherwise.
          )")
      .def("is_isomorphic",
          &mogli::Molecule::is_isomorphic,
          "other"_a,
          R"(
          Test if this molecular graph is isomorphic to another molecular graph.

          Args:
              other (Molecule): Other molecular graph.

          Returns:
              bool. True, if they are isomorphic, false otherwise.
          )")
      .def("print_dot",
          py::overload_cast<>(&mogli::Molecule::print_dot, py::const_),
          R"(
          Export molecular graph to dot (graphviz) format.

          Returns:
              str. Graph in dot format.
          )")
      .def("print_dot",
          py::overload_cast<const mogli::StringVector&>(&mogli::Molecule::print_dot, py::const_),
          "properties"_a,
          R"(
          Export molecular graph to dot (graphviz) format with selected atom properties.

          Args:
              properties (StringVector): Atom properties to print.

          Returns:
              str. Graph in dot format.
          )")
      .def("read_lgf",
          py::overload_cast<const std::string &>(&mogli::Molecule::read_lgf),
          "in"_a,
           R"(
          Import default formatted LGF.

          Args:
              in (str): Molecular graph in LGF format.
          )")
      .def("read_lgf",
          py::overload_cast<const std::string &, const mogli::LGFIOConfig&>(&mogli::Molecule::read_lgf),
          "in"_a, "config"_a,
          R"(
          Import LGF in given format.

          Args:
              in (str):             Molecular graph in LGF format.
              config (LGFIOConfig): LGF formatter.
          )")
      .def("set_property",
          &mogli::Molecule::set_property,
          "node"_a, "property"_a, "value"_a,
          R"(
          Set atom property.

          Args:
              node (Node):                          Atom.
              property (str):                       Property name.
              value (Union[bool, int, float, str]): Property value.
          )")
      .def("split",
           [](mogli::Molecule &self, int max_size, int shell) {
             mogli::SharedPtrVector<mogli::Molecule>::type components;
             self.split(max_size, shell, components);
             return components;
           },
          "max_size"_a, "shell"_a,
          py::return_value_policy::move,
          R"(
          Balanced split of the molecular graph.

          Args:
              max_size (int): Maximal size of the resulting components.
              shell (int):    Shell size (Overlap at the split-regions).

          Returns:
              MoleculeVector. Resulting components.

          Tries to split the molecule into overlapping smaller components of roughly equal size (less or equal max_size),
          )")
      .def("write_lgf",
          py::overload_cast<>(&mogli::Molecule::write_lgf),
          R"(
          Export default formatted LGF.

          Returns:
              str. Molecular graph in default LGF format.
          )")
      .def("write_lgf",
          py::overload_cast<const mogli::LGFIOConfig&>(&mogli::Molecule::write_lgf),
          "config"_a,
          R"(
          Export LGF in given format.

          Args:
              config (LGFIOConfig): LGF formatter.

          Returns:
              str. Molecular graph in given LGF format.
          )");

  py::class_<mogli::Fragment, std::shared_ptr<mogli::Fragment>> fragment(m, "Fragment", molecule);

  fragment.doc() = "Molecular fragment.";

  fragment.def("get_core_atom_count",
          &mogli::Fragment::get_core_atom_count,
          R"(
          Returns the number of core atoms.

          Returns:
              int. Number of core atoms.
          )")
      .def("get_shell_size",
          &mogli::Fragment::get_shell_size,
          R"(
          Returns the shell size of this fragment.

          Returns:
              int. Shell size.
          )")
      .def("is_core",
          &mogli::Fragment::is_core,
          "node"_a,
          R"(
          Test if the given atom is a core atom.

          Args:
              node (Node): Atom.

          Returns:
              bool. True, if the atom is a core atom, false otherwise.
          )")
      .def("print_dot",
          py::overload_cast<>(&mogli::Fragment::print_dot,  py::const_),
          R"(
          Export molecular fragment to dot (graphviz) format.

          Returns:
              str. Fragment in dot format.
          )");

  py::class_<mogli::PeriodicTable> table(m, "PeriodicTable");

  table.doc() = "Periodic table of elements.";

  table.def(
          py::init<>(),
          R"(
          Empty constructor.
          )")
      .def(
          "add",
          &mogli::PeriodicTable::add,
          "num"_a, "name"_a, "color"_a,
          py::return_value_policy::reference_internal,
          R"(
          Add a new element.

          Args:
              num (int):   Element number.
              name (str):  Element type.
              color (str): Color string for dot (graphviz) export.

          Returns:
              PeriodicTable. Update periodic table.
          )")
      .def("add_uncolored",
          &mogli::PeriodicTable::add_uncolored,
           "num"_a, "name"_a,
          py::return_value_policy::reference_internal,
           R"(
          Add a new element.

          Args:
              num (int):   Element number.
              name (str):  Element type.

          Returns:
              PeriodicTable. Update periodic table.
          )");

  py::class_<mogli::Node>(m, "Node")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self);

  // methods

  m.def("are_subgraph_isomorphic",
      [](mogli::Molecule &mol_small, mogli::Molecule &mol_large) {
        mogli::IntToIntMap isomorphism_map;
        bool sub = mogli::are_subgraph_isomorphic(mol_small, mol_large, isomorphism_map);
        return py::make_tuple(sub, isomorphism_map);
      },
      "mol_small"_a, "mol_large"_a,
      py::return_value_policy::move,
      R"(
      Test for subgraph isomorphism of two molecular graphs.

      Args:
          mol_small (Molecule): Smaller molecular graph.
          mol_large (Molecule): Larger molecular graph.

      Returns:
          (bool, IntToIntMap):
               * subgraph (bool):               True, if the smaller graph is a subgraph of the larger graph, false otherwise.
               * isomorphism_map (IntToIntMap): Mapping from the smaller graph to the larger graph.
      )");

  m.def("atomic_fragments",
      [](mogli::Molecule &mol, int shell) {
        mogli::FragmentVector fragments;
        mogli::MatchVector matches;
        mogli::atomic_fragments(mol, fragments, matches, shell);
        return py::make_tuple(fragments, matches);
      },
      "mol"_a, "shell"_a,
      py::return_value_policy::move,
      R"(
      Iterates all atoms of a molecular graph and returns them as fragments a singe-atom core.

      Args:
          mol (Molecule): Molecular graph.
          shell (int):    Shell size of the fragments. Maximal number of bonds from the center atom.

      Returns:
          (FragmentVector, MatchVector)
              * fragments (FragmentVector): Atomic fragments.
              * matches (MatchVector):      Match objects mapping from fragments to the molecular graph.
      )");

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
      "mol1"_a, "mol2"_a, "shell"_a, "min_core_size"_a, "max_core_size"_a,
      "prod_gen"_a, "maximum"_a, "timeout_seconds"_a,
      py::return_value_policy::move,
      R"(
      Computes maximal common fragments of two molecular graphs.

      Args:
          mol1 (Molecule):           First molecular graph.
          mol2 (Molecule):           Second molecular graph.
          shell (int):               Shell size. Maximal number of bonds from any core atom in the fragments.
          min_core_size (int):       Minimal number of core atoms for each fragment.
          max_core_size (int):       Maximal number of core atoms for each fragment.
          prod_gen (GenerationType): Product graph data reduction rule.
          maximum (bool):            If true, reports only the largest fragments.
          timeout_seconds (int):     Timeout in seconds.

      Returns:
          (FragmentVector, MatchVector, MatchVector)
              * fragments (FragmentVector): Maximal common fragments.
              * matches_mol1 (MatchVector): Match objects mapping from fragments to the first molecular graph.
              * matches_mol2 (MatchVector): Match objects mapping from fragments to the second molecular graph.

      The heart and soul of this library. See `this paper <https://doi.org/10.7287/peerj.preprints.3250v1>`_
      for more information. The data reduction rule with the most speedup is GenerationType.UNCON_DEG_1. It is
      recommended to always use this rule, the other rules are mainly for evaluation.
      )");

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
        "mol1"_a, "mol2"_a, "shell"_a, "max_core_size"_a,
        "prod_gen"_a, "maximum"_a, "timeout_seconds"_a,
      py::return_value_policy::move,
        R"(
      Computes maximal common fragments of two molecular graphs.

      Args:
          mol1 (Molecule):           First molecular graph.
          mol2 (Molecule):           Second molecular graph.
          shell (int):               Shell size. Maximal number of bonds from any core atom in the fragments.
          min_core_size (int):       Minimal number of core atoms for each fragment.
          max_core_size (int):       Maximal number of core atoms for each fragment.
          prod_gen (GenerationType): Product graph data reduction rule.
          maximum (bool):            If true, reports only the largest fragments.
          timeout_seconds (int):     Timeout in seconds.

      Returns:
          (FragmentVector, MatchVector, MatchVector)
              * fragments (FragmentVector): Maximal common fragments.
              * matches_mol1 (MatchVector): Match objects mapping from fragments to the first molecular graph.
              * matches_mol2 (MatchVector): Match objects mapping from fragments to the second molecular graph.

      The heart and soul of this library. See `this paper <https://doi.org/10.7287/peerj.preprints.3250v1>`_
      for more information. The data reduction rule with the most speedup is GenerationType.UNCON_DEG_1. It is
      recommended to always use this rule, the other rules are mainly for evaluation.
      )");

  // hashing methods

  m.def("hash_canonization",
      [](const mogli::Canonization & canonization) {
        return py::bytes(mogli::hash_canonization(canonization));
      },
      "canonization"_a,
      py::return_value_policy::move,
      R"(
      Hashes a canonical representation of a molecular graph.

      Args:
          canonization (Canonization): Canonical representation.

      Returns:
          str. Hash.
      )");

  m.def("hash_fcanonization",
      [](const mogli::FragmentCanonization & fcanonization) {
        return py::bytes(mogli::hash_fcanonization(fcanonization));
      },
      "fcanonization"_a,
      py::return_value_policy::move,
      R"(
      Hashes a canonical representation of a molecular fragment.

      Args:
          fcanonization (FragmentCanonization): Canonical representation.

      Returns:
          str. Hash.
      )");

  // serializing methods

  m.def("pack_canonization",
      py::overload_cast<const mogli::Canonization&>(&serialize<const mogli::Canonization&>),
      "obj"_a,
      py::return_value_policy::move,
      R"(
      Serialize a canonical representation of a molecular graph.

      Args:
          obj (Canonization): Canonical representation.

      Returns:
          str. Serialized object.
      )");

  m.def("pack_fcanonization",
      py::overload_cast<const mogli::FragmentCanonization&>(&serialize<const mogli::FragmentCanonization&>),
      "obj"_a,
      py::return_value_policy::move,
      R"(
      Serialize a canonical representation of a molecular fragment.

      Args:
          obj (FragmentCanonization): Canonical representation.

      Returns:
          str. Serialized object.
      )");

  m.def("pack_fragment",
      py::overload_cast<const mogli::Fragment&>(&serialize<mogli::Fragment>),
      "obj"_a,
      py::return_value_policy::move,
      R"(
      Serialize a molecular fragment.

      Args:
          obj (Fragment): Fragment.

      Returns:
          str. Serialized object.
      )");

  m.def("pack_match",
      py::overload_cast<const mogli::Match&>(&serialize<mogli::Match>),
      "obj"_a,
      py::return_value_policy::move,
        R"(
      Serialize a match object.

      Args:
          obj (Match): Match object.

      Returns:
          str. Serialized object.
      )");

  m.def("pack_molecule",
      py::overload_cast<const mogli::Molecule&>(&serialize<mogli::Molecule>),
      "obj"_a,
      py::return_value_policy::move,
      R"(
      Serialize a molecular graph.

      Args:
          obj (Molecule): Molecular graph.

      Returns:
          str. Serialized object.
      )");

  // de-serializing methods

  m.def("unpack_canonization",
      &deserialize<mogli::Canonization>,
      "str"_a,
      py::return_value_policy::move,
      R"(
      Deserialize a canonical representation of a molecular graph.

      Args:
          str (str): Serialized object.

      Returns:
          Canonization. Canonical representation.
      )");

  m.def("unpack_fcanonization",
      &deserialize<mogli::FragmentCanonization>,
      "str"_a,
      py::return_value_policy::move,
      R"(
      Deserialize a canonical representation of a molecular fragment.

      Args:
          str (str): Serialized object.

      Returns:
          FragmentCanonization. Canonical representation.
      )");

  m.def("unpack_fragment",
      &deserialize_ptr<mogli::Fragment>,
      "str"_a,
      py::return_value_policy::reference_internal,
      R"(
      Deserialize a fragment.

      Args:
          str (str): Serialized object.

      Returns:
          Fragment. Molecular fragment.
      )");

  m.def("unpack_match",
      &deserialize<mogli::Match>,
      "str"_a,
      py::return_value_policy::move,
      R"(
      Deserialize a match object.

      Args:
          str (str): Serialized object.

      Returns:
          Match. Match object.
      )");

  m.def("unpack_molecule",
      &deserialize_ptr<mogli::Molecule>,
      "str"_a,
      py::return_value_policy::reference_internal,
      R"(
      Deserialize a molecular graph.

      Args:
          str (str): Serialized object.

      Returns:
          Molecule. Molecular graph.
      )");

}