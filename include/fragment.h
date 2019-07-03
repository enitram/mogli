//
// Created by M. Engler on 06/12/16.
//

#ifndef MOGLI_FRAGMENT_H
#define MOGLI_FRAGMENT_H


#include <deque>
#include <cstring>
#include "product.h"

namespace mogli {

  typedef typename lemon::FilterNodes<const Graph, const NodeToBoolMap> FilterNodes;

  class Fragment : public Molecule {

  private:
    
    typedef std::deque<Node> NodeDeque;
    
    NodeToBoolMap _is_core;
    int _shell_size;
    int _core_node_count;

  public:

    Fragment() :
        Molecule(),
        _is_core(_g, false),
        _shell_size(0),
        _core_node_count(0) {}

    Fragment(const Product &product, const NodeVector &clique, IntToIntMap &g_to_mol1, IntToIntMap &g_to_mol2);

    int get_core_atom_count() const {
      return _core_node_count;
    }

    bool is_core(const Node node) const {
      return _is_core[node];
    }

    void set_core(const Node &u, bool core);

    int get_shell_size() const {
      return _shell_size;
    }

    void set_shell_size(int shell_size) {
      _shell_size = shell_size;
    }

    const std::string print_dot() const override;

    const void print_dot(std::ostream &out) const override;

  private:
    
    void bfs_shell(const Molecule &mol, const Node &v, NodeToNodeMap &shell_g_to_mol, NodeToNodeMap &shell_mol_to_g,
                   NodeToBoolMap &core_nodes, IntSet &shell_ids, NodeToIntMap &shell_min_depth,
                   const Canonization &canon1, const Canonization &canon2, IntToIntMap &g_to_mol1, IntToIntMap &g_to_mol2);

  };

}

#endif //MOGLI_FRAGMENT_H
