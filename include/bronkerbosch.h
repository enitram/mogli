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

#include <chrono>
#include <list>
#include <boost/dynamic_bitset.hpp>
#include "product.h"

namespace mogli {

  class BronKerbosch {
  public:

  private:
    typedef typename boost::dynamic_bitset<> BitSet;
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
    BronKerbosch(const Product& product, unsigned int min_core_size, unsigned int max_core_size, bool maximum);

    void run(int seconds);

    void run(std::chrono::high_resolution_clock::time_point start, long microseconds);

    const NodeVectorVector& getMaxCliques() const {
      return _cliques;
    }

    size_t computeDegeneracy();

  private:

    size_t computeDegeneracy(NodeVector& order);

    void bkPivot(BitSet P, const BitSet & D, const BitSet & R, BitSet X, const BitSet & S,
                 std::chrono::high_resolution_clock::time_point start, long microseconds);

    void report(const BitSet& R);

  };

}

#endif //MOGLI_BRONKERBOSCH_H
