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

    typedef typename Graph::template NodeMap<NodeVector> NodeToNodeVectorMap;
    typedef std::deque<Node> NodeDeque;

    const Molecule &_mol1;
    const Molecule &_mol2;
    NodeToBoolMap _is_shell;
    NodeToBoolMap _is_core;
    NodeToNodeVectorMap _g_to_mol1;
    NodeToNodeVectorMap _g_to_mol2;
    FilterNodes _shell;
    FilterNodes _core;
    int _shell_size;

  public:

    Fragment(const Product &product, const NodeVector &clique) :
        Molecule(),
        _mol1(product.get_mol1()),
        _mol2(product.get_mol2()),
        _is_shell(_g, false),
        _is_core(_g, false),
        _g_to_mol1(_g, NodeVector()),
        _g_to_mol2(_g, NodeVector()),
        _shell(_g, _is_shell),
        _core(_g, _is_core),
        _shell_size(product.get_shell()) {

      for (NodeVector::const_iterator it = clique.begin(), end = clique.end(); it != end; ++it) {
        Node u = product.get_mol1_node(*it);
        Node uv = add_atom(_mol1.get_color(u));
        _g_to_mol1[uv].push_back(u);
        _g_to_mol2[uv].push_back(product.get_mol2_node(*it));
        _is_core[uv] = true;

        const NodePairVector &reductions = product.get_reductions(*it);
        for (NodePairVector::const_iterator it2 = reductions.begin(), end2 = reductions.end(); it2 != end2; ++it2) {
          Node _uv = add_atom(_mol1.get_color(it2->first));
          _g_to_mol1[_uv].push_back(it2->first);
          _g_to_mol2[_uv].push_back(it2->second);
          _is_core[_uv] = true;
        }
      }

      for (NodeIt u(_g); u != lemon::INVALID; ++u) {
        bfs_shell(_g_to_mol1[u].front());
      }

      lemon::ArcLookUp<Graph> arcLookUp(_mol1.get_graph());

      for (NodeIt u(_g); u != lemon::INVALID; ++u) {
        for (NodeIt v(_g); v != lemon::INVALID; ++v) {
          if (arcLookUp(_g_to_mol1[u].front(), _g_to_mol1[v].front()) != lemon::INVALID) {
            _g.addEdge(u,v);
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

    const NodeVector &get_core_mol1_nodes(const Node &node) const {
      return _g_to_mol1[node];
    }

    const NodeVector &get_core_mol2_nodes(const Node &node) const {
      return _g_to_mol2[node];
    }

    void merge(const Fragment &other) {
      // FIXME this can't work, need subgraph isomorphism map!
      for (FilteredNodeIt uv = other.get_core_node_iter(); uv != lemon::INVALID; ++uv) {
        NodeVector other_mol1 = other.get_core_mol1_nodes(uv);
        NodeVector other_mol2 = other.get_core_mol2_nodes(uv);
        _g_to_mol1[uv].insert(_g_to_mol1[uv].end(), other_mol1.begin(), other_mol1.end());
        _g_to_mol2[uv].insert(_g_to_mol2[uv].end(), other_mol2.begin(), other_mol2.end());
      }
    }

  private:

    void bfs_shell(const Node &v) {
      NodeToBoolMap visited(_mol1.get_graph(), false);
      NodeToIntMap depth(_mol1.get_graph(), 0);
      NodeDeque queue;

      queue.push_back(v);
      visited[v] = true;
      while (queue.size() > 0) {
        Node &current = queue.front();
        if (!_is_core[current]) {
          Node uv = add_atom(_mol1.get_color(current));
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
