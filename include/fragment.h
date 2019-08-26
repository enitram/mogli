////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
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

#ifndef MOGLI_FRAGMENT_H
#define MOGLI_FRAGMENT_H


#include <deque>

#include "product.h"

namespace mogli {

  typedef typename lemon::FilterNodes<const Graph, const NodeToBoolMap> FilterNodes;

  /**
   * Molecular fragment.
   */
  class Fragment : public Molecule {

  private:
    
    typedef std::deque<Node> NodeDeque;
    
    NodeToBoolMap _is_core;
    int _shell_size;
    int _core_node_count;

  public:

    /**
     * Empty constructor.
     */
    Fragment() :
        Molecule(),
        _is_core(_g, false),
        _shell_size(0),
        _core_node_count(0) {}

    /**
     * Constructor with custom periodic table.
     *
     * @param[in] periodic_table    Periodic table.
     */
    explicit Fragment(PeriodicTable & periodic_table) :
        Molecule(periodic_table),
        _is_core(_g, false),
        _shell_size(0),
        _core_node_count(0) {}

    /**
     * Constructs a molecular fragment from a clique in the product graph.
     *
     * @param[in] product   Product graph of two molecules.
     * @param[in] clique    Clique in the product graph.
     * @param[in] g_to_mol1 Mapping of product graph nodes to the first molecule.
     * @param[in] g_to_mol2 Mapping of product graph nodes to the second molecule.
     */
    Fragment(const Product &product, const NodeVector &clique, IntToIntMap &g_to_mol1, IntToIntMap &g_to_mol2);

    /**
     * Returns the number of core atoms.
     *
     * @return  Number of core atoms.
     */
    int get_core_atom_count() const {
      return _core_node_count;
    }

    /**
     * Test if the given atom is a core atom.
     *
     * @param[in] node  Atom.
     * @return          True, if the atom is a core atom, false otherwise.
     */
    bool is_core(const Node node) const {
      return _is_core[node];
    }

    /**
     * Set the core atom property for the given atom.
     *
     * @param[in] u     Atom.
     * @param[in] core  Core atom property.
     */
    void set_core(const Node &u, bool core);

    /**
     * Returns the shell size of this fragment.
     *
     * @return  Shell size.
     */
    int get_shell_size() const {
      return _shell_size;
    }

    /**
     * Sets the shell size of this fragment.
     *
     * @param[in] shell_size    Shell size.
     */
    void set_shell_size(int shell_size) {
      _shell_size = shell_size;
    }

    /**
     * Export molecular fragment to dot (graphviz) format.
     *
     * @return  Fragment in dot format.
     */
    const std::string print_dot() const override;

    /**
     * Export molecular fragment to dot (graphviz) format.
     *
     * @param[out] out  Output stream.
     */
    const void print_dot(std::ostream &out) const override;

  private:
    
    void bfs_shell(const Molecule &mol, const Node &v, NodeToNodeMap &shell_g_to_mol, NodeToNodeMap &shell_mol_to_g,
                   NodeToBoolMap &core_nodes, IntSet &shell_ids, NodeToIntMap &shell_min_depth,
                   const Canonization &canon1, const Canonization &canon2, IntToIntMap &g_to_mol1, IntToIntMap &g_to_mol2);

  };

}

#endif //MOGLI_FRAGMENT_H
