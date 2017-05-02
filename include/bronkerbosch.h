/*
 * bronkerbosch.h
 *
 *  Created on: 21-jan-2014
 *      Author: M. El-Kebir
 *
 *  Modified by M. Engler on 02/11/16.
 */

#ifndef MOGLI_BRONKERBOSCH_H
#define MOGLI_BRONKERBOSCH_H

#include "product.h"
#include "boost/graph/vf2_sub_graph_iso.hpp"

namespace mogli {

  class BronKerbosch {
  public:

  private:
    typedef typename Graph::template NodeMap<size_t> NodeToBitMap;
    typedef typename Graph::template NodeMap<BitSet> NodeToBitSetMap;
    typedef std::list<Node> NodeList;
    typedef std::vector<NodeList> NodeListVector;

    const Graph& _g;
    const size_t _n;
    NodeVectorVector _cliques;
    NodeVector _bitToNode;
    NodeToBitMap _nodeToBit;
    NodeToBitSetMap _bitNeighborhood;
    NodeToBitSetMap _restrictedBitNeighborhood;

    unsigned int _min_core_size;
    unsigned int _max_core_size;

    bool _maximum;
    int _current_max;

    const Product &_product;

  public:
    BronKerbosch(const Product& product, unsigned int min_core_size, unsigned int max_core_size, bool maximum)
        : _g(product.get_graph())
        , _product(product)
        , _n(static_cast<size_t>(lemon::countNodes(_g)))
        , _cliques()
        , _bitToNode()
        , _nodeToBit(_g, std::numeric_limits<size_t>::max())
        , _bitNeighborhood(_g, BitSet(_n))
        , _restrictedBitNeighborhood(_g, BitSet(_n))
        , _min_core_size(min_core_size)
        , _max_core_size(max_core_size)
        , _maximum(maximum)
        , _current_max(0) {
      // initialize mappings
      _bitToNode.reserve(_n);
      size_t i = 0;
      for (NodeIt v(_g); v != lemon::INVALID; ++v, ++i) {
        _bitToNode.push_back(v);
        _nodeToBit[v] = i;
      }

      // initialize neighborhoods
      for (NodeIt v(_g); v != lemon::INVALID; ++v, ++i) {
        BitSet& neighborhood = _bitNeighborhood[v];
        for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e) {
          Node w = _g.oppositeNode(v, e);
          neighborhood[_nodeToBit[w]] = 1;
        }
      }

      // initialize restricted neighborhood mapping
      for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
        if (product.is_connectivity_edge(e)) {
          Node u = _g.u(e);
          Node v = _g.v(e);
          _restrictedBitNeighborhood[u][_nodeToBit[v]] = 1;
          _restrictedBitNeighborhood[v][_nodeToBit[u]] = 1;
        }
      }
    }

    void run();

    const NodeVectorVector& getMaxCliques() const {
      return _cliques;
    }

  private:

    size_t computeDegeneracy(NodeVector& order);

    void bkPivot(BitSet P, BitSet D, BitSet R, BitSet X, BitSet S);

    void report(const BitSet& R);

  };

}

#endif //MOGLI_BRONKERBOSCH_H
