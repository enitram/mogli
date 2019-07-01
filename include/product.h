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
#include "orbits.h"
#include <set>
#include <vector>
#include <deque>

namespace mogli {

  typedef typename Graph::template NodeMap<Node> NodeToNodeMap;
  typedef std::vector<std::pair<Node, Node> > NodePairVector;

  typedef std::vector<NodeVector> NodeVectorVector;

  class Product {

  public:

    enum GenerationType {
      NO_OPT = 0,
      UNCON = 1,
      DEG_1 = 2,
      UNCON_DEG_1 = 3,
      SUB = 4,
      UNCON_SUB = 5
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

    typedef typename Graph::template NodeMap<IntSet> NodeToIntSetMap;

    const Molecule& _mol1;
    const Molecule& _mol2;
    const int _shell;

    const GenerationType _gen_type;

    Graph _g;
    NodeToIntMap _node_sizes;
    NodeToNodeMap _g_to_mol1;
    NodeToNodeMap _g_to_mol2;
    EdgeToBoolMap _connectivity;
    NodeToNodePairVectorMap _reductions;
    NodeToCanonizationMap _g_to_mol1_canons;
    NodeToCanonizationMap _g_to_mol2_canons;
    NodeToIntMap _comp_map;
    int _comp_count;
    IntVector _comp_sizes;
    bool _is_complete;

  public:
    Product(const Product& parent, int component)
      : _mol1(parent._mol1)
      , _mol2(parent._mol2)
      , _shell(parent._shell)
      , _gen_type(parent._gen_type)
      , _g()
      , _node_sizes(_g)
      , _reductions(_g)
      , _g_to_mol1(_g)
      , _g_to_mol2(_g)
      , _g_to_mol1_canons(_g)
      , _g_to_mol2_canons(_g)
      , _connectivity(_g)
      , _comp_map(_g)
      , _comp_count(0)
      , _comp_sizes()
      , _is_complete(false) {

      NodeToNodeMap nodes(_g);
      int V = 0;
      for (NodeIt orig(parent._g); orig != lemon::INVALID; ++orig) {
        if (parent._comp_map[orig] == component) {
          Node copy = _g.addNode();
          ++V;
          nodes[copy] = orig;
          _node_sizes[copy] = parent._node_sizes[orig];
          _reductions[copy] = parent._reductions[orig];
          _g_to_mol1[copy] = parent._g_to_mol1[orig];
          _g_to_mol2[copy] = parent._g_to_mol2[orig];
          _g_to_mol1_canons[copy] = parent._g_to_mol1_canons[orig];
          _g_to_mol2_canons[copy] = parent._g_to_mol2_canons[orig];
        }
      }

      lemon::ArcLookUp<Graph> arcLookUp(parent._g);
      int E = 0;
      for (NodeIt v(_g); v != lemon::INVALID; ++v) {
        for (NodeIt u = v; u != lemon::INVALID; ++u) {
          if (v == u)
            continue;

          Edge e = arcLookUp(nodes[u], nodes[v]);
          if (e != lemon::INVALID) {
            _connectivity[_g.addEdge(u,v)] = parent._connectivity[e];
            ++E;
          }
        }
      }
      _is_complete = E == ((V*(V-1))/2);

    }

    Product(const Molecule& mol1, const Molecule& mol2, int shell, GenerationType gen,
            unsigned int min_core_size, unsigned int max_core_size)
      : _mol1(mol1)
      , _mol2(mol2)
      , _shell(shell)
      , _gen_type(gen)
      , _g()
      , _node_sizes(_g)
      , _reductions(_g)
      , _g_to_mol1(_g)
      , _g_to_mol2(_g)
      , _g_to_mol1_canons(_g)
      , _g_to_mol2_canons(_g)
      , _connectivity(_g)
      , _comp_map(_g)
      , _comp_count(1)
      , _is_complete(false) {
      int V, E;
      if ((gen == DEG_1 || gen == UNCON_DEG_1) && shell > 0) {
        V = generate_nodes_deg1();
      } else if ((gen == SUB || gen == UNCON_SUB) && shell > 0) {
        V = generate_nodes_sub();
      } else {
        V = generate_nodes();
      }
      if (gen == UNCON || gen == UNCON_DEG_1 || gen == UNCON_SUB) {
        E = generate_edges_connected(min_core_size, max_core_size);
      } else {
        E = generate_edges();
      }
      _is_complete = E == ((V*(V-1))/2);

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

    const int get_components() const {
      return _comp_count;
    }

    const int get_component_size(int component) const {
      return _comp_sizes[component];
    }

    const int get_shell() const {
      return _shell;
    }

    const bool is_complete() const {
      return _is_complete;
    }

    const GenerationType get_gen_type() const {
      return _gen_type;
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

    const Canonization& get_mol1_canon(const Node& uv) const {
      return _g_to_mol1_canons[uv];
    }

    const Canonization& get_mol2_canon(const Node& uv) const {
      return _g_to_mol2_canons[uv];
    }

    const int get_clique_size(const NodeVector& clique) const {
      int size = 0;
      for (auto& v : clique) {
        size += _node_sizes[v];
      }
      return size;
    }

    bool is_connectivity_edge(Edge e) const {
      return _connectivity[e];
    }

    const std::string print_dot() const {
      std::stringstream buffer;
      print_dot(buffer);
      return buffer.str();
    }

    const std::string print_dot(const StringVector &properties) const {
      std::stringstream buffer;
      print_dot(buffer, properties);
      return buffer.str();
    }

    const void print_dot(std::ostream& out) const {
      // header
      out << "graph G {" << std::endl
          << "\toverlap=scale" << std::endl;

      // nodes
      for (NodeIt uv(_g); uv != lemon::INVALID; ++uv) {
        Node u = _g_to_mol1[uv];
        Node v = _g_to_mol2[uv];

        out << "\t" << _g.id(uv);
        out << "[style=\"filled\",fillcolor=" << _mol1.get_color_name(u);
        out << ",label=\"" << _mol1.get_id(u) << "x" << _mol2.get_id(v) << "\"]";
        out << std::endl;
      }

      // edges
      for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
        out << "\t" << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
        if (_connectivity[e]) {
          out <<  std::endl;
        } else {
          out << "[style=\"dashed\"]" << std::endl;
        }
      }

      out << "}" << std::endl;
    }

    const void print_dot(std::ostream &out, const StringVector &properties) const {
      // header
      out << "graph G {" << std::endl
          << "\toverlap=scale" << std::endl;

      // nodes
      for (NodeIt uv(_g); uv != lemon::INVALID; ++uv) {
        Node u = _g_to_mol1[uv];
        Node v = _g_to_mol2[uv];

        out << "\t" << _g.id(uv);
        out << "[style=\"filled\",fillcolor=" << _mol1.get_color_name(u);
        out << ",label=\"";
        bool first = true;
        for (std::string prop : properties) {
          if (!first) {
            out << ",";
          } else {
            first = false;
          }
          Any value1 = _mol1.get_property(u, prop);
          if (std::holds_alternative<bool>(value1)) {
            out << std::get<bool>(value1);
          } else if (std::holds_alternative<int>(value1)) {
            out << std::get<int>(value1);
          } else if (std::holds_alternative<double>(value1)) {
            out << std::get<double>(value1);
          } else if (std::holds_alternative<std::string>(value1)) {
            out << std::get<std::string>(value1);
          }
          out << "x";
          Any value2 = _mol2.get_property(v, prop);
          if (std::holds_alternative<bool>(value2)) {
            out << std::get<bool>(value2);
          } else if (std::holds_alternative<int>(value2)) {
            out << std::get<int>(value2);
          } else if (std::holds_alternative<double>(value2)) {
            out << std::get<double>(value2);
          } else if (std::holds_alternative<std::string>(value2)) {
            out << std::get<std::string>(value2);
          }
        }
        out << "\"]";
        out << std::endl;
      }

      // edges
      for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
        out << "\t" << _g.id(_g.u(e)) << " -- " << _g.id(_g.v(e));
        if (_connectivity[e]) {
          out <<  std::endl;
        } else {
          out << "[style=\"dashed\"]" << std::endl;
        }
      }

      out << "}" << std::endl;
    }
    
  private:

    Node add_node(const Node& u, const Node& v);

    int generate_nodes();

    int generate_nodes_deg1();

    int generate_nodes_sub();

    int generate_edges();

    int generate_edges_connected(unsigned int min_core_size, unsigned int max_core_size);

    void determine_degrees(const Graph& g, IntToNodeMap& deg_to_node, NodeToIntMap& deg);

    void generate_subgraph_canonization(const Molecule &mol, const Node &v, NodeToCanonizationMap &map);

    void generate_subgraph(const Molecule &mol, const Node &v, NodeToIntSetMap &neighborhoods, IntToNodeMap &sizes);

    void bfs(const Molecule &mol, const Node &v, NodeToBoolMap &filter);

    void bfs_neighbors(const Molecule &mol, const Node &v, IntSet &neighbors, int& size);

    void bfs_subgraph(const Molecule &mol, const Node &product_node, const Node &root_node, const IntSet &root_neighbors,
                      const NodeToIntSetMap &neighborhoods, const NodeVector &order1, const NodeVector &order2,
                      ShortToNodeVectorPairMap &current_reductions);

  };

} // namespace mogli

#endif // MOGLI_PRODUCT_H
