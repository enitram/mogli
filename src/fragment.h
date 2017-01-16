//
// Created by M. Engler on 06/12/16.
//

#ifndef MOGLI_FRAGMENT_H
#define MOGLI_FRAGMENT_H


#include <deque>
#include "product.h"

namespace mogli {

  typedef typename lemon::FilterNodes<const Graph, const NodeToBoolMap> FilterNodes;
  typedef typename FilterNodes::NodeIt FilteredNodeIt;
  typedef typename FilterNodes::EdgeIt FilteredEdgeIt;
  typedef typename FilterNodes::IncEdgeIt FilteredIncEdgeIt;

  class Fragment : public Molecule {

  private:

    typedef typename Graph::template NodeMap<StringVector> NodeToStringVectorMap;
    typedef std::deque<Node> NodeDeque;

    const Molecule& _mol1;
    const Molecule& _mol2;
    NodeToBoolMap _is_shell;
    NodeToBoolMap _is_core;
    NodeToStringVectorMap _g_to_mol1;
    NodeToStringVectorMap _g_to_mol2;
    FilterNodes _shell;
    FilterNodes _core;
    int _shell_size;
    std::string _unp;

  public:

    Fragment(const Molecule& mol1, const Molecule& mol2, const std::string unique_node_property) :
        Molecule(),
        _mol1(mol1),
        _mol2(mol2),
        _g_to_mol1(_g),
        _g_to_mol2(_g),
        _is_shell(_g, false),
        _is_core(_g, false),
        _shell(_g, _is_shell),
        _core(_g, _is_core),
        _shell_size(0),
        _unp(unique_node_property) {}

    Fragment(const Product &product, const NodeVector &clique, const std::string unique_node_property) :
        Molecule(),
        _mol1(product.get_mol1()),
        _mol2(product.get_mol2()),
        _is_shell(_g, false),
        _is_core(_g, false),
        _g_to_mol1(_g, StringVector()),
        _g_to_mol2(_g, StringVector()),
        _shell(_g, _is_shell),
        _core(_g, _is_core),
        _shell_size(product.get_shell()),
        _unp(unique_node_property) {

      NodeToNodeMap shell_nodes(_g);
      NodeToBoolMap core_nodes(_mol1.get_graph(), false);

      for (NodeVector::const_iterator it = clique.begin(), end = clique.end(); it != end; ++it) {
        Node u = product.get_mol1_node(*it);
        Node uv = add_atom(_mol1.get_color(u));
        _g_to_mol1[uv].push_back(_mol1.get_string_property(u, _unp));
        _g_to_mol2[uv].push_back(_mol2.get_string_property(product.get_mol2_node(*it), _unp));
        _is_core[uv] = true;
        core_nodes[u] = true;

        const NodePairVector &reductions = product.get_reductions(*it);
        for (NodePairVector::const_iterator it2 = reductions.begin(), end2 = reductions.end(); it2 != end2; ++it2) {
          Node _uv = add_atom(_mol1.get_color(it2->first));
          _g_to_mol1[_uv].push_back(_mol1.get_string_property(it2->first, _unp));
          _g_to_mol2[_uv].push_back(_mol1.get_string_property(it2->second, _unp));
          _is_core[_uv] = true;
        }
      }

      for (NodeIt u(_g); u != lemon::INVALID; ++u) {
        bfs_shell(_mol1.get_node_by_string_property(_unp, _g_to_mol1[u].front()), shell_nodes, core_nodes);
      }

      lemon::ArcLookUp<Graph> arcLookUp(_mol1.get_graph());

      for (NodeIt uv1(_g); uv1 != lemon::INVALID; ++uv1) {
        for (NodeIt uv2 = uv1; uv2 != lemon::INVALID; ++uv2) {
          if (uv1 == uv2)
            continue;
          Node u1 = _is_core[uv1] ? _mol1.get_node_by_string_property(_unp, _g_to_mol1[uv1].front())
                                  : shell_nodes[uv1];
          Node u2 = _is_core[uv2] ? _mol1.get_node_by_string_property(_unp, _g_to_mol1[uv2].front())
                                  : shell_nodes[uv2];
          if (arcLookUp(u1, u2) != lemon::INVALID) {
            _g.addEdge(uv1, uv2);
          }
        }
      }

    }

    const FilteredNodeIt get_core_node_iter() const {
      return FilteredNodeIt(_core);
    }

    const FilteredNodeIt get_shell_node_iter() const {
      return FilteredNodeIt(_shell);
    }

    const FilteredEdgeIt get_core_edge_iter() const {
      return FilteredEdgeIt(_core);
    }

    const FilteredEdgeIt get_shell_edge_iter() const {
      return FilteredEdgeIt(_shell);
    }

    const FilteredIncEdgeIt get_core_inc_edge_iter(const Node &node) const {
      return FilteredIncEdgeIt(_core, node);
    }

    const FilteredIncEdgeIt get_shell_inc_edge_iter(const Node &node) const {
      return FilteredIncEdgeIt(_shell, node);
    }

    const Molecule& get_mol1() const {
      return _mol1;
    }

    const Molecule &get_mol2() const {
      return _mol2;
    }

    void get_core_mol1_unp(const Node &node, StringVector &vector) const {
      for (StringVector::const_iterator it = _g_to_mol1[node].begin(), end = _g_to_mol1[node].end(); it != end; ++it) {
        vector.push_back(*it);
      }
    }

    void get_core_mol2_unp(const Node &node, StringVector &vector) const {
      for (StringVector::const_iterator it = _g_to_mol2[node].begin(), end = _g_to_mol2[node].end(); it != end; ++it) {
        vector.push_back(*it);
      }
    }

    const Node get_core_mol1_node(const Node &node) const {
      return _mol1.get_node_by_string_property(_unp, _g_to_mol1[node].front());
    }

    const Node get_core_mol2_node(const Node &node) const {
      return _mol2.get_node_by_string_property(_unp, _g_to_mol2[node].front());
    }

    void add_core_mol1_unp(const Node &node, std::string unique_node_property) {
      _g_to_mol1[node].push_back(unique_node_property);
    }

    void add_core_mol2_unp(const Node &node, std::string unique_node_property) {
      _g_to_mol2[node].push_back(unique_node_property);
    }

    void merge(const Fragment &other) {
      // TODO merging fragments of different molecules -> ptr_list of molecules (then we can also make an empty constructor)?
      // FIXME this can't work, need subgraph isomorphism map!
      for (FilteredNodeIt uv = other.get_core_node_iter(); uv != lemon::INVALID; ++uv) {
        StringVector other_mol1_nodes, other_mol2_nodes;
        other.get_core_mol1_unp(uv, other_mol1_nodes);
        other.get_core_mol2_unp(uv, other_mol2_nodes);
//        for (NodeVector::const_iterator it = other_mol1_nodes.begin(), end = other_mol1_nodes.end(); it != end; ++it) {
//          _g_to_mol1[uv].push_back(_mol1.get_string_property(*it, _unp));
//        }
//        for (NodeVector::const_iterator it = other_mol2_nodes.begin(), end = other_mol2_nodes.end(); it != end; ++it) {
//          _g_to_mol2[uv].push_back(_mol2.get_string_property(*it, _unp));
//        }

        _g_to_mol1[uv].insert(_g_to_mol1[uv].end(), other_mol1_nodes.begin(), other_mol1_nodes.end());
        _g_to_mol2[uv].insert(_g_to_mol2[uv].end(), other_mol2_nodes.begin(), other_mol2_nodes.end());
      }
    }

    bool is_shell(const Node &u) const {
      return _is_shell[u];
    }

    void set_shell(const Node &u, const bool shell) {
      _is_shell[u] = shell;
      _is_core[u] = !shell;
    }

    int get_shell_size() const {
      return _shell_size;
    }

    void set_shell_size(int shell_size) {
      _shell_size = shell_size;
    }

  private:

    void bfs_shell(const Node &v, NodeToNodeMap &shell_nodes, NodeToBoolMap &core_nodes) {
      NodeToIntMap depth(_mol1.get_graph(), 0);
      NodeToBoolMap visited(_mol1.get_graph(), false);
      NodeDeque queue;

      queue.push_back(v);
      visited[v] = true;
      while (queue.size() > 0) {
        Node &current = queue.front();
        if (!core_nodes[current]) {
          Node uv = add_atom(_mol1.get_color(current));
          shell_nodes[uv] = current;
          _is_shell[uv] = true;
        }

        if (depth[current] < _shell_size) {
          for (IncEdgeIt e = _mol1.get_inc_edge_iter(current); e != lemon::INVALID; ++e) {
            Node w = _mol1.get_opposite_node(current, e);
            if (!visited[w]) {
              visited[w] = true;
              depth[w] = depth[current]+1;
              queue.push_back(w);
            }
          }
        }
        queue.pop_front();
      }
    }

  };

}

#endif //MOGLI_FRAGMENT_H
