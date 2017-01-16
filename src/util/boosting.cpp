//
// Created by martin on 24/10/16.
//

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "canonization.h"
#include "packing.h"
#include "isomorphism.h"
#include "mcf.h"

using namespace boost::python;
using namespace mogli;

inline boost::python::object pass_through(boost::python::object const& o) { return o; }

template<class Klass, class KlassIter>
struct iterator_wrappers {

  static Klass
  next(KlassIter& o) {
    Klass next = o;
    if (next == lemon::INVALID) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      throw_error_already_set();
    }
    ++o;
    return next;
  }

  static void
  wrap(const char* python_name) {
    class_<KlassIter>(python_name)
        .def("next", next)
        .def("__iter__", pass_through);
  }

};

BOOST_PYTHON_MODULE(libmogli) {

  def("pack_canonization", pack_canonization);
  def("unpack_canonization", unpack_canonization);
  def("hash_canonization", hash_canonization);
  def("pack_molecule", pack_molecule);
  def("unpack_molecule", unpack_molecule);
  def("pack_fragment", pack_fragment);
  def("unpack_fragment", unpack_fragment);

  class_<ShortVector>("ShortVector")
      .def("__iter__", iterator<ShortVector>())
      .def(vector_indexing_suite<ShortVector>());

  class_<LongVector>("LongVector")
      .def("__iter__", iterator<LongVector>())
      .def(vector_indexing_suite<LongVector>());

  class_<StringVector>("StringVector")
      .def("__iter__", iterator<StringVector>())
      .def(vector_indexing_suite<StringVector>());

  class_<NodeVector>("NodeVector")
      .def("__iter__", iterator<NodeVector>())
      .def(vector_indexing_suite<NodeVector>());

  class_<FragmentVector, boost::noncopyable>("FragmentVector")
      .def("__iter__", iterator<FragmentVector, return_internal_reference<> >())
      .def("__len__", &FragmentVector::size);

  class_<Canonization>("Canonization", init<const Molecule&>())
      .def("get_colors", &Canonization::get_colors, return_value_policy<return_by_value>())
      .def("get_canonization", &Canonization::get_canonization, return_value_policy<return_by_value>())
      .def("get_node_order", &Canonization::get_node_order, return_value_policy<return_by_value>());

  enum_<Product::GenerationType>("GenerationType")
      .value("NO_OPT", Product::GenerationType::NO_OPT)
      .value("DEG_1", Product::GenerationType::DEG_1)
      .value("SUB", Product::GenerationType::SUB);

  const Node (Molecule::*add_atom1)(std::string) = &Molecule::add_atom;
  const Node (Molecule::*add_atom2)(unsigned short) = &Molecule::add_atom;

  void (Molecule::*set_property1)(Node, std::string, bool) = &Molecule::set_property;
  void (Molecule::*set_property2)(Node, std::string, int) = &Molecule::set_property;
  void (Molecule::*set_property3)(Node, std::string, double) = &Molecule::set_property;
  void (Molecule::*set_property4)(Node, std::string, char*) = &Molecule::set_property;
  void (Molecule::*set_property5)(Node, std::string, std::string) = &Molecule::set_property;

  class_<Molecule, boost::noncopyable>("Molecule")
      .def("add_atom", add_atom1)
      .def("add_atom", add_atom2)
      .def("add_edge", &Molecule::add_edge)
      .def("get_node_iter", &Molecule::get_node_iter)
      .def("get_edge_iter", &Molecule::get_edge_iter)
      .def("get_inc_edge_iter", &Molecule::get_inc_edge_iter)
      .def("get_opposite_node", &Molecule::get_opposite_node)
      .def("get_u", &Molecule::get_u)
      .def("get_v", &Molecule::get_v)
      .def("get_atom_count", &Molecule::get_atom_count)
      .def("get_node_by_id", &Molecule::get_node_by_id)
      .def("get_id", &Molecule::get_id)
      .def("get_color", &Molecule::get_color)
      .def("get_iacm_element", &Molecule::get_iacm_element)
      .def("get_chem_element", &Molecule::get_chem_element)
      .def("is_connected", &Molecule::is_connected)
      .def("read_lgf", &Molecule::read_lgf)
      .def("add_bool_property", &Molecule::add_bool_property)
      .def("add_int_property", &Molecule::add_int_property)
      .def("add_double_property", &Molecule::add_double_property)
      .def("add_string_property", &Molecule::add_string_property)
      .def("get_bool_properties", &Molecule::get_bool_properties)
      .def("get_int_properties", &Molecule::get_int_properties)
      .def("get_double_properties", &Molecule::get_double_properties)
      .def("get_string_properties", &Molecule::get_string_properties)
      .def("set_property", set_property1)
      .def("set_property", set_property2)
      .def("set_property", set_property3)
      .def("set_property", set_property4)
      .def("set_property", set_property5)
      .def("get_bool_property", &Molecule::get_bool_property)
      .def("get_int_property", &Molecule::get_int_property)
      .def("get_double_property", &Molecule::get_double_property)
      .def("get_string_property", &Molecule::get_string_property)
      .def("get_node_by_string_property", &Molecule::get_node_by_string_property);

  class_<Fragment, bases<Molecule>, boost::noncopyable>("Fragment", init<const Molecule&, const Molecule&, const std::string>())
      .def("get_core_node_iter", &Fragment::get_core_node_iter)
      .def("get_shell_node_iter", &Fragment::get_shell_node_iter)
      .def("get_core_edge_iter", &Fragment::get_core_edge_iter)
      .def("get_shell_edge_iter", &Fragment::get_shell_edge_iter)
      .def("get_core_inc_edge_iter", &Fragment::get_core_inc_edge_iter)
      .def("get_shell_inc_edge_iter", &Fragment::get_shell_inc_edge_iter)
//      .def("get_core_mol1_unp", &Fragment::get_core_mol1_unp)
//      .def("get_core_mol2_unp", &Fragment::get_core_mol2_unp)
      .def("get_core_mol1_node", &Fragment::get_core_mol1_node)
      .def("get_core_mol2_node", &Fragment::get_core_mol2_node);

  bool (*isomorph1)(Molecule&, Molecule&) = &are_isomorphic;
  bool (*isomorph2)(Canonization&, Canonization&) = &are_isomorphic;

  def("are_isomorphic", isomorph1);
  def("are_isomorphic", isomorph2);

  def("maximal_common_fragments", maximal_common_fragments);

  class_<Node>("Node")
      .def(self == self)
      .def(self != self)
      .def(self < self);

  class_<Edge>("Edge")
      .def(self == self)
      .def(self != self)
      .def(self < self);

  iterator_wrappers<const Node, NodeIt>().wrap("NodeIt");
  iterator_wrappers<const Edge, EdgeIt>().wrap("EdgeIt");
  iterator_wrappers<const Edge, IncEdgeIt>().wrap("IncEdgeIt");

  iterator_wrappers<const Node, FilteredNodeIt>().wrap("FilteredNodeIt");
  iterator_wrappers<const Edge, FilteredEdgeIt>().wrap("FilteredEdgeIt");
  iterator_wrappers<const Edge, FilteredIncEdgeIt>().wrap("FilteredIncEdgeIt");

}
