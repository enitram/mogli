//
// Created by martin on 10/21/16.
//

#ifndef MOGLI_PACKING_H
#define MOGLI_PACKING_H

#include <msgpack.hpp>
#include <sstream>
#include "../fcanonization.h"
#include "../canonization.h"
#include "../fragment.h"
#include "../match.h"

namespace msgpack {
  MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
    namespace adaptor {

      using namespace mogli;

      template<>
      struct pack<Canonization> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, Canonization const& v) const {
          return o.pack_array(3)
              .pack(v.get_colors())
              .pack(v.get_canonization())
              .pack(v.get_node_order());
        }
      };

      template<>
      struct convert<Canonization> {
        msgpack::object const& operator()(msgpack::object const& o, Canonization& v) const {
          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 3) throw msgpack::type_error();
          v = Canonization(
              o.via.array.ptr[0].as<ShortVector>(),
              o.via.array.ptr[1].as<LongVector>(),
              o.via.array.ptr[2].as<ShortVector>()
          );
          return o;
        }
      };

      template<>
      struct pack<FragmentCanonization> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, FragmentCanonization const& v) const {
          return o.pack_array(4)
              .pack(v.get_colors())
              .pack(v.get_canonization())
              .pack(v.get_node_order())
              .pack(v.get_core_nodes());
        }
      };

      template<>
      struct convert<FragmentCanonization> {
        msgpack::object const& operator()(msgpack::object const& o, FragmentCanonization& v) const {
          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 4) throw msgpack::type_error();
          v = FragmentCanonization(
              o.via.array.ptr[0].as<ShortVector>(),
              o.via.array.ptr[1].as<LongVector>(),
              o.via.array.ptr[2].as<ShortVector>(),
              o.via.array.ptr[3].as<BoolVector>()
          );
          return o;
        }
      };

      template<>
      struct pack<boost::any> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, boost::any const& v) const {
          if (v.type() == typeid(bool)) {
            o.pack_array(2).pack(0).pack(boost::any_cast<bool>(v));
          } else if (v.type() == typeid(int)) {
            o.pack_array(2).pack(1).pack(boost::any_cast<int>(v));
          } else if (v.type() == typeid(double)) {
            o.pack_array(2).pack(2).pack(boost::any_cast<double>(v));
          } else if (v.type() == typeid(std::string)) {
            o.pack_array(2).pack(3).pack(boost::any_cast<std::string>(v));
          } else {
            throw msgpack::type_error();
          }
          return o;
        }
      };

      template<>
      struct convert<boost::any> {

        msgpack::object const& operator()(msgpack::object const& o, boost::any& v) const {
          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 2) throw msgpack::type_error();

          msgpack::object *arr = o.via.array.ptr;
          int type = arr[0].as<int>();
          switch (type)  {
            case 0: v = boost::any(arr[1].as<bool>());
              break;
            case 1: v = boost::any(arr[1].as<int>());
              break;
            case 2: v = boost::any(arr[1].as<double>());
              break;
            case 3: v = boost::any(arr[1].as<std::string>());
              break;
            default:
              throw msgpack::type_error();
          }
          return o;
        }
      };

      template<>
      struct pack<Molecule> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, Molecule const& v) const {
          const Graph &g = v.get_graph();
          const int n = v.get_atom_count();

          o.pack_array(4);
          o.pack_array(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            o.pack_array(2).pack(v.get_id(u)).pack(v.get_color(u));
          }

          std::vector<std::pair<int, int> > edges;
          for (EdgeIt e = EdgeIt(g); e != lemon::INVALID; ++e) {
            edges.push_back(std::make_pair(v.get_id(g.v(e)), v.get_id(g.u(e))));
          }
          o.pack(edges);


          StringVector properties;
          v.get_properties(properties);

          o.pack(properties);
          o.pack_array(properties.size());
          for (std::string& property : properties) {
            o.pack_array(n);
            for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
              o.pack_array(2).pack(v.get_id(u)).pack(v.get_property(u, property));
            }
          }
          return o;
        }
      };

      msgpack::object const& convert_molecule(msgpack::object const& o, Molecule& mol) {
        if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
        if (o.via.array.size != 4) throw msgpack::type_error();

        msgpack::object *arr = o.via.array.ptr;
        for (int i = 0; i < 4; ++i)
          if (arr[i].type != msgpack::type::ARRAY) throw msgpack::type_error();

        std::vector<std::pair<int, unsigned short> > atoms = arr[0].as<std::vector<std::pair<int, unsigned short> > >();

        for (std::vector<std::pair<int, unsigned short> >::const_iterator it = atoms.begin(), end = atoms.end(); it != end; ++it) {
          mol.add_atom(it->first, it->second);
        }

        std::vector<std::pair<int, int> > edges = arr[1].as<std::vector<std::pair<int, int> > >();
        for (std::vector<std::pair<int, int> >::const_iterator it = edges.begin(), end = edges.end(); it != end; ++it) {
          Node u = mol.get_node_by_id(it->first);
          Node v = mol.get_node_by_id(it->second);
          mol.add_edge(u, v);
        }

        StringVector props = arr[2].as<StringVector>();

        std::vector<std::vector<std::pair<int, boost::any> > > properties = arr[3].as<std::vector<std::vector<std::pair<int, boost::any> > > >();
        if (props.size() != properties.size()) throw msgpack::type_error();
        for (int i = 0; i < props.size(); ++i) {
          for (auto& pair : properties[i]) {
            mol.set_property(mol.get_node_by_id(pair.first), props[i], pair.second);
          }
        }

        return o;
      }

      template<>
      struct convert<Molecule> {

        msgpack::object const& operator()(msgpack::object const& o, Molecule& v) const {
          return convert_molecule(o, v);
        }
      };

      template<>
      struct pack<Fragment> {
        template <typename Stream>
        packer<Stream>& operator()(msgpack::packer<Stream>& o, Fragment const& v) const {
          const Molecule &mol = v;
          o.pack_array(3).pack(mol);

          const Graph& g = v.get_graph();
          int n = v.get_atom_count();

          o.pack(v.get_shell_size());
          o.pack_array(n);
          for (NodeIt u = NodeIt(g); u != lemon::INVALID; ++u) {
            o.pack_array(2).pack(v.get_id(u)).pack(v.is_core(u));
          }

          return o;
        }
      };

      template<>
      struct convert<Fragment> {
        msgpack::object const &operator()(msgpack::object const &o, Fragment &v) const {

          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 3) throw msgpack::type_error();
          msgpack::object *arr = o.via.array.ptr;

          convert_molecule(arr[0], v);

          v.set_shell_size(arr[1].as<int>());

          std::vector<std::pair<int, bool> > shell = arr[2].as<std::vector<std::pair<int, bool> > >();
          for (std::vector<std::pair<int, bool> >::const_iterator it = shell.begin(), end = shell.end(); it != end; ++it) {
            v.set_core(v.get_node_by_id(it->first), it->second);
          }

          return o;
        }
      };

      template<>
      struct pack<Match> {
        template<typename Stream>
        packer<Stream> &operator()(msgpack::packer<Stream> &o, Match const &v) const {
          o.pack_array(2).pack(v.get_frag_to_mol()).pack(v.get_merged_frag_to_mol());
          return o;
        }
      };

      template<>
      struct convert<Match> {
        msgpack::object const &operator()(msgpack::object const &o, Match &v) const {

          if (o.type != msgpack::type::ARRAY) throw msgpack::type_error();
          if (o.via.array.size != 2) throw msgpack::type_error();

          msgpack::object *arr = o.via.array.ptr;

          IntToIntMap ftm = arr[0].as<IntToIntMap>();
          IntToIntMapVector mftm = arr[1].as<IntToIntMapVector>();

          for (IntToIntMap::const_iterator it = ftm.begin(), end = ftm.end(); it != end; ++it) {
            v.add_frag_to_mol(it->first, it->second);
          }
          for (IntToIntMapVector::iterator it = mftm.begin(), end = mftm.end(); it != end; ++it) {
            v.add_merged_frag_to_mol(*it);
          }

          return o;
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

  std::string pack_fcanonization(const FragmentCanonization &obj) {
    std::stringstream buffer;
    pack(buffer, obj);
    return buffer.str();
  }

  void unpack_fcanonization(std::string str, FragmentCanonization &obj) {
    unpacked msg_unpacked = unpack(str.data(), str.size());
    object msg_object = msg_unpacked.get();
    msg_object.convert(obj);
  }

  std::string hash_fcanonization(const FragmentCanonization &obj) {
    std::stringstream buffer;
    packer<std::stringstream> p(buffer);
    p.pack_array(3).pack(obj.get_colors()).pack(obj.get_canonization()).pack(obj.get_core_nodes());
    return buffer.str();
  }

  std::string pack_molecule(const Molecule &obj) {
    std::stringstream buffer;
    pack(buffer, obj);
    return buffer.str();
  }

  std::string pack_molecule(const boost::shared_ptr<Molecule> &obj) {
    return pack_molecule(*obj);
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

  std::string pack_fragment(const boost::shared_ptr<Fragment> &obj) {
    return pack_fragment(*obj);
  }

  void unpack_fragment(std::string str, Fragment &obj) {
    unpacked msg_unpacked = unpack(str.data(), str.size());
    object msg_object = msg_unpacked.get();
    msg_object.convert(obj);
  }

  std::string pack_match(const Match &obj) {
    std::stringstream buffer;
    pack(buffer, obj);
    return buffer.str();
  }

  void unpack_match(std::string str, Match &obj) {
    unpacked msg_unpacked = unpack(str.data(), str.size());
    object msg_object = msg_unpacked.get();
    msg_object.convert(obj);
  }

}

#endif //MOGLI_PACKING_H
