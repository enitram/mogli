////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
//    Copyright (C) 2014       M. El-Kebir                                                                            //
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

#ifndef MOGLI_PRODUCT_H
#define MOGLI_PRODUCT_H

#include <deque>

#include "molecule.h"
#include "canonization.h"


namespace mogli {

  typedef typename Graph::template NodeMap<Node> NodeToNodeMap;
  typedef std::vector<std::pair<Node, Node> > NodePairVector;

  typedef std::vector<NodeVector> NodeVectorVector;

  class Product {

  public:

    /**
     * @brief Product graph data reduction rules.
     *
     * The data reduction rule with the most speedup is UNCON_DEG_1. It is recommended to always use this rule,
     * the other rules are mainly for evaluation.
     */
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
    Product(const Product& parent, int component);

    Product(const Molecule& mol1, const Molecule& mol2, int shell, GenerationType gen,
            unsigned int min_core_size);

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

    std::string print_dot() const;

    std::string print_dot(const StringVector &properties) const;

    void print_dot(std::ostream& out) const;

    void print_dot(std::ostream &out, const StringVector &properties) const;
    
  private:

    Node add_node(const Node& u, const Node& v);

    int generate_nodes();

    int generate_nodes_deg1();

    int generate_nodes_sub();

    int generate_edges();

    int generate_edges_connected(unsigned int min_core_size);

    static void determine_degrees(const Graph& g, IntToNodeMap& deg_to_node, NodeToIntMap& deg);

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
