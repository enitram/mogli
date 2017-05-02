//
// Created by M. Engler on 02/11/16.
//

#include "../include/bronkerbosch.h"

using namespace mogli;

void BronKerbosch::run() {

  NodeVector order;
  computeDegeneracy(order);

  BitSet mask(_n);

  for (typename NodeVector::const_iterator it = order.begin(); it != order.end(); ++it) {
    Node v = *it;
    const BitSet& N_v = _bitNeighborhood[v];
    const BitSet& Nc_v = _restrictedBitNeighborhood[v];
    const BitSet Nd_v = N_v - Nc_v;

    // ~mask includes v but we're fine as N_v, Nc_v and Nd_v exclude v
    BitSet P = Nc_v & ~mask;
    BitSet D = Nd_v & ~mask;
    BitSet X = Nc_v & mask;
    BitSet S = Nd_v & mask;

    BitSet R(_n);
    R.set(_nodeToBit[v]);

    bkPivot(P, D, R, X, S);
    mask.set(_nodeToBit[v]);
  }

}

size_t BronKerbosch::computeDegeneracy(NodeVector& order)
{
  // Requires O(|V| + |E|) time
  order.clear();

  typedef typename Graph::template NodeMap<size_t> DegNodeMap;
  typedef typename Graph::template NodeMap<typename NodeList::iterator> NodeListItMap;

  NodeToBoolMap present(_g, true);
  DegNodeMap deg(_g, 0);
  size_t maxDeg = 0;
  NodeListItMap it(_g);

  // compute node degrees, O(|E|) time
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    size_t d = 0;
    for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e, ++d);
    deg[v] = d;
    if (d > maxDeg) maxDeg = d;
  }

  // fill T, O(d) time
  NodeListVector T(maxDeg + 1, NodeList());
  for (NodeIt v(_g); v != lemon::INVALID; ++v) {
    size_t d = deg[v];
    T[d].push_front(v);
    it[v] = T[d].begin();
  }

  size_t degeneracy = 0;

  // O(|V|) time, Eppstein et al. (2010)
  const size_t n = T.size();
  size_t i = 0;
  while (i < n) {
    NodeList& l = T[i];
    if (T[i].size() > 0) {
      Node v = l.front();
      l.pop_front();
      order.push_back(v);
      present[v] = false;
      if (deg[v] > degeneracy) {
        degeneracy = deg[v];
      }

      for (IncEdgeIt e(_g, v); e != lemon::INVALID; ++e) {
        Node w = _g.oppositeNode(v, e);
        if (present[w]) {
          size_t deg_w = deg[w];
          typename NodeList::iterator it_w = it[w];

          T[deg_w - 1].splice(T[deg_w - 1].begin(), T[deg_w], it_w);
          deg[w]--;
        }
      }

      i = 0;
    } else {
      ++i;
    }
  }

  return degeneracy;
}

void BronKerbosch::bkPivot(BitSet P, BitSet D,
                           BitSet R,
                           BitSet X, BitSet S) {
  // all sets are pairwise disjoint
  assert((P & D).none());
  assert((P & R).none());
  assert((P & X).none());
  assert((P & S).none());
  assert((D & R).none());
  assert((D & X).none());
  assert((D & S).none());
  assert((R & X).none());
  assert((R & S).none());
  assert((X & S).none());

  // Reports maximal c-cliques in P \cup R (but not in X and S)
  BitSet P_cup_X = P | X;
  if (P_cup_X.none()) {
    report(R);
  } else {
    // choose a pivot u from (P | X) s.t |P & N(u)| is maximum, Tomita et al. (2006)
    size_t maxBitCount = 0;
    Node max_u = lemon::INVALID;
    for (size_t i = 0; i < P.size(); ++i) {
      if (P_cup_X[i]) {
        Node u = _bitToNode[i];
        BitSet P_cap_Nu = P & _bitNeighborhood[u];
        size_t s = P_cap_Nu.count();
        if (s >= maxBitCount) {
          max_u = u;
          maxBitCount = s;
        }
      }
    }

    assert(max_u != lemon::INVALID);
    BitSet P_diff_Nu = P - _bitNeighborhood[max_u];
    for (size_t i = 0; i < P.size(); ++i) {
      if (P_diff_Nu[i]) {
        Node v = _bitToNode[i];
        const BitSet& N_v = _bitNeighborhood[v];
        const BitSet& Nc_v = _restrictedBitNeighborhood[v];

        BitSet P_ = P | ( D & Nc_v );
        BitSet D_ = D - Nc_v;

        BitSet R_ = R;
        R_[_nodeToBit[v]] = 1;

        BitSet X_ = X | ( S & Nc_v );
        BitSet S_ = S - Nc_v;

        // report all maximal cliques in ( (P | N[v]) & R) \ (X & N[v]) )
        bkPivot(P_ & N_v,
                D_ & N_v,
                R_,
                X_ & N_v,
                S_ & N_v);
        P[i] = 0;
        X[i] = 1;
      }
    }
  }

  /* Invariants:
   * - R is a c-clique
   * - Each node v in P is adjacent to all nodes in R
   *   and c-adjacent to a node in R
   *   (nodes in P are used to extend R, upon usage it's moved to X)
   * - Each node v in D is adjacent to all nodes in R
   *   but not c-adjacent to any node in R
   * - Each node v in X is adjacent to all nodes in R
   *   (maximal c-cliques containing R \cup v have already been reported)
   */
}

void BronKerbosch::report(const BitSet& R) {
  NodeVector clique;
  for (BitSet::size_type it = R.find_first(); it != BitSet::npos; it = R.find_next(it)) {
    clique.push_back(_bitToNode[static_cast<int>(it)]);
  }
  int size = _product.get_clique_size(clique);
  if (_min_core_size <= size && size <= _max_core_size) {
    if (!_maximum) {
      _cliques.push_back(clique);
    } else if (size > _current_max) {
      _cliques.clear();
      _current_max = size;
      _cliques.push_back(clique);
    } else if (size == _current_max) {
      _cliques.push_back(clique);
    }
  }
}

