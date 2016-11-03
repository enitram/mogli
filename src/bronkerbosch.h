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
#include <boost/dynamic_bitset.hpp>

namespace mogli {

  typedef std::vector<NodeVector> NodeMatrix;

  class BronKerbosch {
  public:

    typedef enum {
      BK_CLASSIC,
      BK_PIVOT,
      BK_PIVOT_DEGENERACY
    } SolverType;

  private:
    typedef typename Graph::template NodeMap<size_t> NodeToBitMap;
    typedef boost::dynamic_bitset<> BitSet;
    typedef typename Graph::template NodeMap<BitSet> NodeToBitSetMap;

    const Graph& _g;
    const size_t _n;
    NodeMatrix _cliques;
    NodeVector _bitToNode;
    NodeToBitMap _nodeToBit;
    NodeToBitSetMap _bitNeighborhood;

  public:
    BronKerbosch(const Graph& g)
        : _g(g)
        , _n(static_cast<size_t>(lemon::countNodes(_g)))
        , _cliques()
        , _bitToNode()
        , _nodeToBit(g, std::numeric_limits<size_t>::max())
        , _bitNeighborhood(g, BitSet(_n)) {
      // initialize mappings
      _bitToNode.reserve(_n);
      size_t i = 0;
      for (NodeIt v(_g); v != lemon::INVALID; ++v, ++i)
      {
        _bitToNode.push_back(v);
        _nodeToBit[v] = i;
      }

      // initialize neighborhoods
      for (NodeIt v(_g); v != lemon::INVALID; ++v, ++i)
      {
        BitSet& neighborhood = _bitNeighborhood[v];
        for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e)
        {
          Node w = _g.oppositeNode(v, e);
          neighborhood[_nodeToBit[w]] = 1;
        }
      }
    }

    virtual void run(SolverType type);

//    void print(std::ostream& out) const;

    size_t getNumberOfMaxCliques() const {
      return _cliques.size();
    }

    const NodeMatrix& getMaxCliques() const {
      return _cliques;
    }

  private:

//    size_t computeDegeneracy(NodeVector& order);

    void report(const BitSet& R);

//    void printBitSet(const BitSet& S, std::ostream& out) const;

    /// Classic Bron-Kerbosch algorithm without pivoting
    ///
    /// Reports maximal cliques in P \cup R (but not in X)
    void bkClassic(BitSet P, BitSet R, BitSet X);
//    void bkPivot(BitSet P, BitSet R, BitSet X);
//    void bkDegeneracy(const NodeVector& order);
  };

}

#endif //MOGLI_BRONKERBOSCH_H
