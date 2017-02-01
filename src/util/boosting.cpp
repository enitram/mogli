//
// Created by martin on 24/10/16.
//

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "canonization.h"
#include "packing.h"
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

  std::string (*pack_fragment1)(const Fragment&) = &pack_fragment;
  std::string (*pack_fragment2)(const boost::shared_ptr<Fragment>&) = &pack_fragment;

  def("pack_canonization", pack_canonization);
  def("unpack_canonization", unpack_canonization);
  def("hash_canonization", hash_canonization);
  def("pack_fcanonization", pack_fcanonization);
  def("unpack_fcanonization", unpack_fcanonization);
  def("hash_fcanonization", hash_fcanonization);
  def("pack_molecule", pack_molecule);
  def("unpack_molecule", unpack_molecule);
  def("pack_fragment", pack_fragment1);
  def("pack_fragment", pack_fragment2);
  def("unpack_fragment", unpack_fragment);
  def("pack_match", pack_match);
  def("unpack_match", unpack_match);

  class_<BoolVector>("BoolVector")
      .def("__iter__", iterator<BoolVector>())
      .def(vector_indexing_suite<BoolVector>());

  class_<ShortVector>("ShortVector")
      .def("__iter__", iterator<ShortVector>())
      .def(vector_indexing_suite<ShortVector>());

  class_<IntVector>("IntVector")
      .def("__iter__", iterator<IntVector>())
      .def(vector_indexing_suite<IntVector>());

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

  class_<MatchVector, boost::noncopyable>("MatchVector")
      .def("__iter__", iterator<MatchVector>())
      .def("__len__", &MatchVector::size);

  class_<Canonization>("Canonization", init<const Molecule&>())
      .def("get_colors", &Canonization::get_colors, return_value_policy<return_by_value>())
      .def("get_canonization", &Canonization::get_canonization, return_value_policy<return_by_value>())
      .def("get_node_order", &Canonization::get_node_order, return_value_policy<return_by_value>())
      .def("is_isomorphic", &Canonization::is_isomorphic);

  class_<FragmentCanonization, bases<Canonization>>("FragmentCanonization", init<const Fragment&>())
      .def("get_core_nodes", &FragmentCanonization::get_core_nodes, return_value_policy<return_by_value>())
      .def("is_isomorphic", &FragmentCanonization::is_isomorphic);

  enum_<Product::GenerationType>("GenerationType")
      .value("NO_OPT", Product::GenerationType::NO_OPT)
      .value("DEG_1", Product::GenerationType::DEG_1)
      .value("SUB", Product::GenerationType::SUB);

  const Node (Molecule::*add_atom1)(std::string) = &Molecule::add_atom;
  const Node (Molecule::*add_atom2)(int, std::string) = &Molecule::add_atom;
  const Node (Molecule::*add_atom3)(unsigned short) = &Molecule::add_atom;
  const Node (Molecule::*add_atom4)(int, unsigned short) = &Molecule::add_atom;

  void (Molecule::*set_property1)(Node, std::string, bool) = &Molecule::set_property;
  void (Molecule::*set_property2)(Node, std::string, int) = &Molecule::set_property;
  void (Molecule::*set_property3)(Node, std::string, double) = &Molecule::set_property;
  void (Molecule::*set_property4)(Node, std::string, char*) = &Molecule::set_property;
  void (Molecule::*set_property5)(Node, std::string, std::string) = &Molecule::set_property;

  void (Molecule::*read_lgf1)(const std::string &) = &Molecule::read_lgf;
  void (Molecule::*read_lgf2)(const std::string &, const std::string, const std::string) = &Molecule::read_lgf;

  const std::string (Molecule::*print_dot1)(const StringVector&) const = &Molecule::print_dot;
  const std::string (Fragment::*print_dot2)() const = &Fragment::print_dot;

  class_<Molecule, boost::noncopyable>("Molecule", init<>())
      .def("add_atom", add_atom1)
      .def("add_atom", add_atom2)
      .def("add_atom", add_atom3)
      .def("add_atom", add_atom4)
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
      .def("is_isomorphic", &Molecule::is_isomorphic)
      .def("read_lgf", read_lgf1)
      .def("read_lgf", read_lgf2)
      .def("print_dot", print_dot1)
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
      .def("get_string_property", &Molecule::get_string_property);

  class_<Fragment, boost::shared_ptr<Fragment>, bases<Molecule>, boost::noncopyable>("Fragment", init<>())
      .def("get_core_atom_count", &Fragment::get_core_atom_count)
      .def("is_core", &Fragment::is_core)
      .def("print_dot", print_dot2);

  class_<Match>("Match", init<>())
      .def("frag_to_mol", &Match::frag_to_mol)
      .def("merged_frag_to_mol", &Match::merged_frag_to_mol)
      .def("map_ids", &Match::map_ids)
      .def(self == self)
      .def(self != self);

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

}
