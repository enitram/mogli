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

#ifndef MOGLI_CANONIZATION_H
#define MOGLI_CANONIZATION_H

#include "molecule.h"

namespace mogli {


  /**
   * Canonical representation of a molecular graph.
   */
  class Canonization {

  public:

    /**
     * Initialize an empty canonization.
     */
    Canonization() = default;

    /**
     * Create a canonical representation of a molecular graph.
     *
     * @param[in] mol   Molecular graph.
     */
    explicit Canonization(const Molecule& mol) : _colors(), _canonization(), _node_order() {
      init(mol);
    }

    /**
     * Create a canonical representation of a subgraph of a molecular graph.
     *
     * @param[in] mol       Molecular graph.
     * @param[in] filter    Subgraph filter.
     * @param[in] root      Special atom, always first in canonical representation.
     */
    Canonization(const Molecule& mol, const NodeToBoolMap& filter, const Node& root) :
        _colors(), _canonization(), _node_order() {
      init(mol, filter, root);
    }

    /**
     * Move constructor.
     *
     * @param[in] _colors           Element numbers.
     * @param[in] _canonization     Canonical representations.
     * @param[in] _node_order       Atom IDs.
     */
    Canonization(ShortVector _colors, LongVector _canonization, ShortVector _node_order) :
        _colors(std::move(_colors)),
        _canonization(std::move(_canonization)),
        _node_order(std::move(_node_order)) {}

    /**
     * Returns the element numbers of the atoms in canonical order.
     *
     * @return  Element numbers.
     */
    const ShortVector &get_colors() const {
      return _colors;
    }

    /**
     * Returns the canonical representations of the atoms.
     *
     * @return  Canonical representations.
     */
    const LongVector &get_canonization() const {
      return _canonization;
    }

    /**
     * Returns the atom IDs in canonical order.
     *
     * @return  Atom IDs.
     */
    const ShortVector &get_node_order() const {
      return _node_order;
    }

    /**
     * @brief Isomorphism test.
     *
     * The element matching function determines when two atoms are matching. Usually, two atoms match if they are
     * of the same element (and thus have the same element numbers).
     *
     * @param[in] other     Other canonization.
     * @param[in] matcher   Element number matching function.
     * @return              True, if isomorphic to other canonization, false otherwise.
     */
    const bool is_isomorphic(const Canonization &other, const ElementMatcher & matcher = &default_matcher) const;

  protected:

    typedef std::set<unsigned short> ShortSet;
    typedef typename Graph::template NodeMap<bool> NodeToBoolMap;
    typedef typename Graph::template NodeMap<int> NodeToIntMap;
    typedef typename std::map<unsigned short, NodeVector> ShortToNodeVectorMap;
    typedef typename lemon::FilterNodes<const Graph, const NodeToBoolMap> FilterNodes;
    typedef typename FilterNodes::NodeIt FilteredNodeIt;
    typedef typename FilterNodes::EdgeIt FilteredEdgeIt;
    typedef typename FilterNodes::IncEdgeIt FilteredIncEdgeIt;
    typedef typename FilterNodes::NodeMap<bool> FilteredNodeToBoolMap;

    ShortVector _colors;
    LongVector _canonization;
    ShortVector _node_order;

    void init(const Molecule& mol);

    void init(const Molecule& mol, const NodeToBoolMap& filter, const Node& root);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             NodeToBoolMap& visited, ShortSet& colorSet, ShortToNodeVectorMap& colorMap);

    void dfs(const Node& current, const Node& last, const Molecule& mol,
             const FilterNodes& subgraph, NodeToBoolMap& visited, ShortSet& colorSet,
             ShortToNodeVectorMap& colorMap, unsigned int& node_count);

    void canonNauty(const Molecule& mol,
                    const ShortSet &colorSet,
                    const ShortToNodeVectorMap &colorMap,
                    unsigned int atom_count);

    void canonNauty(const Molecule& mol,
                    const FilterNodes& subgraph,
                    const ShortSet &colorSet,
                    const ShortToNodeVectorMap &colorMap,
                    const Node& root,
                    unsigned int atom_count);

  };

}

#endif //MOGLI_CANONIZATION_H
