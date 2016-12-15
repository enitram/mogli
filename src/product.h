/*
 * product.h
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 *
 *  Modified by M. Engler on 26/10/16.
 */

#ifndef MOGLI_PRODUCT_H
#define MOGLI_PRODUCT_H

#include "molecule.h"
#include "canonization.h"
#include "isomorphism.h"
#include "orbits.h"
#include <set>
#include <vector>
#include <deque>

namespace mogli {

  typedef typename Graph::template NodeMap<Node> NodeToNodeMap;
  typedef std::vector<std::pair<Node, Node> > NodePairVector;
  typedef std::vector<NodePairVector> NodePairVectorVector;

  typedef std::vector<NodeVector> NodeVectorVector;

  class Product {

  public:

    // TODO SUB rule without connection to root node -> check connection after matching to decide if shell or core!
    // TODO SUB rule new implementation (slower, for correctness check only): build all product nodes, check against every other neighborhood

    enum GenerationType {
      NO_OPT = 0,
      DEG_1 = 1,
      SUB = 2
    };

  private:
    typedef typename Graph::template EdgeMap<bool> EdgeToBoolMap;
    typedef typename Graph::template NodeMap<int> NodeToIntMap;

    typedef std::deque<Node> NodeDeque;

    typedef std::pair<NodeVector, NodeVector> NodeVectorPair;
    typedef std::map<unsigned short, NodeVectorPair> ShortToNodeVectorPairMap;

    typedef std::multimap<int, Node> IntToNodeMap;

    typedef typename Graph::template NodeMap<NodePairVector> NodeToNodePairVectorMap;
    typedef typename Graph::template NodeMap<Canonization> NodeToCanonizationMap;

    typedef typename Graph::template NodeMap<BitSet> NodeToBitSetMap;

    typedef std::vector<int> IntVector;
    typedef std::pair<IntVector, IntVector> IntVectorPair;
    typedef std::vector<IntVectorPair> IntVectorPairVector;

    const Molecule& _mol1;
    const Molecule& _mol2;
    const int _shell;

    const GenerationType _gen_type;

    Graph _g;
    NodeToNodeMap _g_to_mol1;
    NodeToNodeMap _g_to_mol2;
    EdgeToBoolMap _connectivity;
    NodeToNodePairVectorMap _reductions;

  public:
    Product(const Molecule& mol1, const Molecule& mol2, int shell, GenerationType gen)
      : _mol1(mol1)
      , _mol2(mol2)
      , _shell(shell)
      , _gen_type(gen)
      , _g()
      , _reductions(_g)
      , _g_to_mol1(_g)
      , _g_to_mol2(_g)
      , _connectivity(_g) {
      NodeVectorVector nodes;
      if (gen == DEG_1 && shell > 0) {
        generate_nodes_deg1();
      } else if (gen == SUB && shell > 0) {
        generate_nodes_sub();
      } else {
        generate_nodes();
      }
      generate_edges();
    }

    const Molecule& get_mol1() const {
      return _mol1;
    }

    const Molecule& get_mol2() const {
      return _mol2;
    }

    const Graph& get_graph() const {
      return _g;
    }

    const int get_shell() const {
      return _shell;
    }

    const Node& get_mol1_node(const Node& uv) const {
      return _g_to_mol1[uv];
    }

    const Node& get_mol2_node(const Node& uv) const {
      return _g_to_mol2[uv];
    }

    const NodePairVector& get_reductions(const Node& uv) const {
      return _reductions[uv];
    }

    bool is_connectivity_edge(Edge e) const {
      return _connectivity[e];
    }

    const void print_dot(std::ostream& out, const StringVector& properties = {}) const {
      // header
      out << "graph G {" << std::endl
          << "\toverlap=scale" << std::endl
          << "\tlayout=neato" << std::endl;

      StringVector bool_props, int_props, double_props, string_props;
      if (properties.size() > 0) {
        StringVector keys1;
        StringVector keys2;
        _mol1.get_bool_properties(keys1);
        _mol2.get_bool_properties(keys2);
        check_properties(properties, keys1, keys2, bool_props);
        _mol1.get_int_properties(keys1);
        _mol2.get_int_properties(keys2);
        check_properties(properties, keys1, keys2, int_props);
        _mol1.get_double_properties(keys1);
        _mol2.get_double_properties(keys2);
        check_properties(properties, keys1, keys2, double_props);
        _mol1.get_string_properties(keys1);
        _mol2.get_string_properties(keys2);
        check_properties(properties, keys1, keys2, string_props);
      }

      // nodes
      for (NodeIt uv(_g); uv != lemon::INVALID; ++uv) {
        Node u = _g_to_mol1[uv];
        Node v = _g_to_mol2[uv];

        out << "\t" << _g.id(uv);
        if (properties.size() > 0) {
          out << "[style=filled,fillcolor=" << _mol1.get_chem_color(u);
          out << ",label=\"";
          bool first = true;
          for (std::vector<std::string>::const_iterator it = string_props.begin(), end = string_props.end(); it != end; ++it) {
            if (first) {
              out << _mol1.get_string_property(u, *it) << "," << _mol2.get_string_property(v, *it);
              first = false;
            } else {
              out << "\\n" << _mol1.get_string_property(u, *it) << "," << _mol2.get_string_property(v, *it);
            }
          }
          for (std::vector<std::string>::const_iterator it = double_props.begin(), end = double_props.end(); it != end; ++it) {
            if (first) {
              out << _mol1.get_double_property(u, *it) << "," << _mol2.get_double_property(v, *it);
              first = false;
            } else {
              out << "\\n" << _mol1.get_double_property(u, *it) << "," << _mol2.get_double_property(v, *it);
            }
          }
          for (std::vector<std::string>::const_iterator it = int_props.begin(), end = int_props.end(); it != end; ++it) {
            if (first) {
              out << _mol1.get_int_property(u, *it) << "," << _mol2.get_int_property(v, *it);
              first = false;
            } else {
              out << "\\n" << _mol1.get_int_property(u, *it) << "," << _mol2.get_int_property(v, *it);
            }
          }
          for (std::vector<std::string>::const_iterator it = bool_props.begin(), end = bool_props.end(); it != end; ++it) {
            if (first) {
              out << _mol1.get_bool_property(u, *it) << "," << _mol2.get_bool_property(v, *it);
              first = false;
            } else {
              out << "\\n" << _mol1.get_bool_property(u, *it) << "," << _mol2.get_bool_property(v, *it);
            }
          }
          out << "\"]";
        }
        out << std::endl;
      }

      // edges
      for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
        out << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
        if (!is_connectivity_edge(e)) {
          out << " [style=dashed]";
        }
        out << std::endl;
      }

      out << "}" << std::endl;
    }

//    void get_node_mapping(const NodeVectorVector &cliques, NodePairVectorVector &mapping);
//
//    void prune_cliques(NodeVectorVector &cliques);
    
  private:

    const void inline check_properties(const StringVector& from, const StringVector& keys1, const StringVector& keys2, StringVector& to) const {
      for (std::vector<std::string>::const_iterator it = from.begin(), end = from.end(); it != end; ++it) {
        if (std::find(keys1.begin(), keys1.end(), *it) != keys1.end() && std::find(keys2.begin(), keys2.end(), *it) != keys2.end()) {
          to.push_back(*it);
        }
      }
    }

    Node add_node(const Node& u, const Node& v);

    void generate_nodes();

    void generate_nodes_deg1();

    void generate_nodes_sub();

    void generate_edges();

    void determine_degrees(const Graph& g, IntToNodeMap& deg_to_node, NodeToIntMap& deg);

    void generate_subgraph_canonization(const Molecule &mol, const Node &v, NodeToCanonizationMap &map);

    void generate_subgraph(const Molecule &mol, const Node &v, NodeToBitSetMap &neighborhoods, IntToNodeMap &sizes);

    void bfs(const Molecule &mol, const Node &v, NodeToBoolMap &filter);

    void bfs_neighbors(const Molecule &mol, const Node &v, BitSet &neighbors, int& size);

    void bfs_subgraph(const Molecule &mol, const Node &product_node, const Node &root_node, const BitSet &root_neighbors,
                      const NodeToBitSetMap &neighborhoods, const NodeVector &order1, const NodeVector &order2,
                      ShortToNodeVectorPairMap &current_reductions);

  };

} // namespace mogli

#endif // MOGLI_PRODUCT_H
