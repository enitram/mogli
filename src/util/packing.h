//
// Created by martin on 10/21/16.
//

#ifndef MOGLI_PACKING_H
#define MOGLI_PACKING_H

#include <msgpack.hpp>
#include <sstream>
#include "canonization.h"
#include "fragment.h"

namespace msgpack {
  MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
    namespace adaptor {

      using namespace mogli;

      template<>
      struct pack<Canonization> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, Canonization const& v) const {
          return o.pack_array(2).pack(v.get_colors()).pack(v.get_canonization()).pack(v.get_node_order());
        }
      };

      template<>
      struct convert<Canonization> {
        msgpack::object const& operator()(msgpack::object const& o, Canonization& v) const {
          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 2) throw msgpack::type_error();
          v = Canonization(
              o.via.array.ptr[0].as<ShortVector>(),
              o.via.array.ptr[1].as<LongVector>(),
              o.via.array.ptr[2].as<ShortVector>()
          );
          return o;
        }
      };

      template<>
      struct pack<Molecule> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, Molecule const& v) const {
          const Graph &g = v.get_graph();
          const int n = v.get_atom_count();

          o.pack_array(10);
          o.pack_array(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            o.pack_array(2).pack(g.id(u)).pack(v.get_color(u));
          }

          o.pack_map(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            o.pack(g.id(u));
            std::vector<int> edges;
            for (IncEdgeIt e = IncEdgeIt(g, u); e != lemon::INVALID; ++e) {
              edges.push_back(g.id(g.oppositeNode(u, e)));
            }
            o.pack(edges);
          }

          StringVector bool_props, int_props, double_props, string_props;
          v.get_bool_properties(bool_props);
          v.get_int_properties(int_props);
          v.get_double_properties(double_props);
          v.get_string_properties(string_props);

          o.pack(bool_props).pack(int_props).pack(double_props).pack(string_props);

          o.pack_array(bool_props.size());
          for (StringVector::const_iterator it = bool_props.begin(), end = bool_props.end(); it != end; ++it) {
            o.pack_array(n);
            for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
              o.pack_array(2).pack(g.id(u)).pack(v.get_bool_property(u, *it));
            }
          }
          o.pack_array(int_props.size());
          for (StringVector::const_iterator it = int_props.begin(), end = int_props.end(); it != end; ++it) {
            o.pack_array(n);
            for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
              o.pack_array(2).pack(g.id(u)).pack(v.get_int_property(u, *it));
            }
          }
          o.pack_array(double_props.size());
          for (StringVector::const_iterator it = double_props.begin(), end = double_props.end(); it != end; ++it) {
            o.pack_array(n);
            for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
              o.pack_array(2).pack(g.id(u)).pack(v.get_double_property(u, *it));
            }
          }
          o.pack_array(string_props.size());
          for (StringVector::const_iterator it = string_props.begin(), end = string_props.end(); it != end; ++it) {
            o.pack_array(n);
            for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
              o.pack_array(2).pack(g.id(u)).pack(v.get_string_property(u, *it));
            }
          }
          return o;
        }
      };

      msgpack::object const& convert_molecule(msgpack::object const& o, Molecule& v, std::map<int, Node> &nodes) {
        if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
        if (o.via.array.size != 10) throw msgpack::type_error();

        msgpack::object *arr = o.via.array.ptr;
        if (arr[0].type != msgpack::type::ARRAY) throw msgpack::type_error();
        if (arr[1].type != msgpack::type::MAP) throw msgpack::type_error();
        for (int i = 2; i < 10; ++i)
          if (arr[i].type != msgpack::type::ARRAY) throw msgpack::type_error();

        std::vector<std::pair<int, unsigned short> > atoms = arr[0].as<std::vector<std::pair<int, unsigned short> > >();

        for (std::vector<std::pair<int, unsigned short> >::const_iterator it = atoms.begin(), end = atoms.end(); it != end; ++it) {
          nodes[it->first] = v.add_atom(it->second);
        }

        std::map<int, std::vector<int> > edges = arr[1].as<std::map<int, std::vector<int> > >();
        for (std::map<int, std::vector<int> >::const_iterator it = edges.begin(), end = edges.end(); it != end; ++it) {
          Node u = nodes[it->first];
          for (std::vector<int>::const_iterator it2 = it->second.begin(), end2 = it->second.end(); it2 != end2; ++it2) {
            v.add_edge(u, nodes[*it2]);
          }
        }

        StringVector bool_props = arr[2].as<StringVector>();
        StringVector int_props = arr[3].as<StringVector>();
        StringVector double_props = arr[4].as<StringVector>();
        StringVector string_props = arr[5].as<StringVector>();

        for (StringVector::const_iterator it = bool_props.begin(), end = bool_props.end(); it != end; ++it) {
          v.add_bool_property(*it);
        }
        for (StringVector::const_iterator it = int_props.begin(), end = int_props.end(); it != end; ++it) {
          v.add_int_property(*it);
        }
        for (StringVector::const_iterator it = double_props.begin(), end = double_props.end(); it != end; ++it) {
          v.add_double_property(*it);
        }
        for (StringVector::const_iterator it = string_props.begin(), end = string_props.end(); it != end; ++it) {
          v.add_string_property(*it);
        }

        std::vector<std::vector<std::pair<int, bool> > > bool_properties = arr[6].as<std::vector<std::vector<std::pair<int, bool> > > >();
        if (bool_props.size() != bool_properties.size()) throw msgpack::type_error();
        for (int i = 0; i < bool_props.size(); ++i) {
          v.add_bool_property(bool_props[i]);
          for (std::vector<std::pair<int, bool> >::const_iterator it = bool_properties[i].begin(), end = bool_properties[i].end(); it != end; ++it) {
            v.set_property(nodes[it->first], bool_props[i], it->second);
          }
        }

        std::vector<std::vector<std::pair<int, int> > > int_properties = arr[7].as<std::vector<std::vector<std::pair<int, int> > > >();
        if (int_props.size() != int_properties.size()) throw msgpack::type_error();
        for (int i = 0; i < int_props.size(); ++i) {
          v.add_int_property(int_props[i]);
          for (std::vector<std::pair<int, int> >::const_iterator it = int_properties[i].begin(), end = int_properties[i].end(); it != end; ++it) {
            v.set_property(nodes[it->first], int_props[i], it->second);
          }
        }

        std::vector<std::vector<std::pair<int, double> > > double_properties = arr[8].as<std::vector<std::vector<std::pair<int, double> > > >();
        if (double_props.size() != double_properties.size()) throw msgpack::type_error();
        for (int i = 0; i < double_props.size(); ++i) {
          v.add_double_property(double_props[i]);
          for (std::vector<std::pair<int, double> >::const_iterator it = double_properties[i].begin(), end = double_properties[i].end(); it != end; ++it) {
            v.set_property(nodes[it->first], double_props[i], it->second);
          }
        }

        std::vector<std::vector<std::pair<int, std::string> > > string_properties = arr[9].as<std::vector<std::vector<std::pair<int, std::string> > > >();
        if (string_props.size() != string_properties.size()) throw msgpack::type_error();
        for (int i = 0; i < string_props.size(); ++i) {
          v.add_string_property(string_props[i]);
          for (std::vector<std::pair<int, std::string> >::const_iterator it = string_properties[i].begin(), end = string_properties[i].end(); it != end; ++it) {
            v.set_property(nodes[it->first], string_props[i], it->second);
          }
        }

        return o;
      }

      template<>
      struct convert<Molecule> {

        msgpack::object const& operator()(msgpack::object const& o, Molecule& v) const {
          std::map<int, Node> nodes;
          return convert_molecule(o, v, nodes);
        }
      };

      template<>
      struct pack<Fragment> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, Fragment const& v) const {
          const Molecule &mol = v;
          o.pack_array(5).pack(mol);

          const Graph& g = v.get_graph();
          int n = v.get_atom_count();

          o.pack(v.get_shell_size());
          o.pack_array(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            o.pack_array(2).pack(g.id(u)).pack(v.is_shell(u));
          }

          const Molecule &mol1 = v.get_mol1();
          o.pack_array(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            StringVector unp;
            v.get_core_mol1_unp(u, unp);
            o.pack_array(2).pack(g.id(u)).pack(unp);
          }

          const Molecule &mol2 = v.get_mol2();
          o.pack_array(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            StringVector unp;
            v.get_core_mol2_unp(u, unp);
            o.pack_array(2).pack(g.id(u)).pack(unp);
          }

          return o;
        }
      };

      template<>
      struct convert<Fragment> {
        msgpack::object const &operator()(msgpack::object const &o, Fragment &v) const {
          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 5) throw msgpack::type_error();
          msgpack::object *arr = o.via.array.ptr;

          Molecule &mol = v;
          std::map<int, Node> nodes;
          convert_molecule(arr[0], mol, nodes);

          v.set_shell_size(arr[1].as<int>());

          std::vector<std::pair<int, bool> > shell = arr[2].as<std::vector<std::pair<int, bool> > >();
          for (std::vector<std::pair<int, bool> >::const_iterator it = shell.begin(), end = shell.end(); it != end; ++it) {
            v.set_shell(nodes[it->first], it->second);
          }

          const Molecule &mol1 = v.get_mol1();
          std::vector<std::pair<int, std::vector<std::string> > > mol1_ids = arr[3].as<std::vector<std::pair<int, std::vector<std::string> > > >();
          for (std::vector<std::pair<int, std::vector<std::string> > >::const_iterator it = mol1_ids.begin(), end = mol1_ids.end(); it != end; ++it) {
            for (std::vector<std::string>::const_iterator it2 = it->second.begin(), end2 = it->second.end(); it2 != end2; ++it2) {
              v.add_core_mol1_unp(nodes[it->first], *it2);
            }
          }

          const Molecule &mol2 = v.get_mol1();
          std::vector<std::pair<int, std::vector<std::string> > > mol2_ids = arr[4].as<std::vector<std::pair<int, std::vector<std::string> > > >();
          for (std::vector<std::pair<int, std::vector<std::string> > >::const_iterator it = mol2_ids.begin(), end = mol2_ids.end(); it != end; ++it) {
            for (std::vector<std::string>::const_iterator it2 = it->second.begin(), end2 = it->second.end(); it2 != end2; ++it2) {
              v.add_core_mol2_unp(nodes[it->first], *it2);
            }
          }

        }
      };

    } // namespace adaptor
  } // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack

namespace mogli {

  using namespace msgpack;

  std::string pack_canonization(const Canonization &obj) {
    std::stringstream buffer;
    pack(buffer, obj);
    return buffer.str();
  }

  void unpack_canonization(std::string str, Canonization &obj) {
    unpacked msg_unpacked = unpack(str.data(), str.size());
    object msg_object = msg_unpacked.get();
    msg_object.convert(obj);
  }

  std::string hash_canonization(const Canonization &obj) {
    std::stringstream buffer;
    packer<std::stringstream> p(buffer);
    p.pack_array(2).pack(obj.get_colors()).pack(obj.get_canonization());
    return buffer.str();
  }

  std::string pack_molecule(const Molecule &obj) {
    std::stringstream buffer;
    pack(buffer, obj);
    return buffer.str();
  }

  void unpack_molecule(std::string str, Molecule &obj) {
    unpacked msg_unpacked = unpack(str.data(), str.size());
    object msg_object = msg_unpacked.get();
    msg_object.convert(obj);
  }

  std::string pack_fragment(const Fragment &obj) {
    std::stringstream buffer;
    pack(buffer, obj);
    return buffer.str();
  }

  void unpack_fragment(std::string str, Fragment &obj) {
    unpacked msg_unpacked = unpack(str.data(), str.size());
    object msg_object = msg_unpacked.get();
    msg_object.convert(obj);
  }

}

#endif //MOGLI_PACKING_H
