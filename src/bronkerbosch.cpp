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

#include "bronkerbosch.h"


using namespace mogli;

BronKerbosch::BronKerbosch(const mogli::Product &product, unsigned int min_core_size, unsigned int max_core_size, bool maximum) :
    _g(product.get_graph())
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
      neighborhood[_nodeToBit[w]] = true;
    }
  }

  // initialize restricted neighborhood mapping
  for (EdgeIt e(_g); e != lemon::INVALID; ++e) {
    if (product.is_connectivity_edge(e)) {
      Node u = _g.u(e);
      Node v = _g.v(e);
      _restrictedBitNeighborhood[u][_nodeToBit[v]] = true;
      _restrictedBitNeighborhood[v][_nodeToBit[u]] = true;
    }
  }
}

void BronKerbosch::run(int seconds) {
  std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
  long microseconds = 1000l * seconds;
  run(start, microseconds);
}

bool BronKerbosch::run(std::chrono::high_resolution_clock::time_point start, long microseconds) {

  // the productgraph is a complete graph with a spanning tree of c-edges, we don't actually need to run BK!
  if ((_product.get_gen_type() == Product::GenerationType::UNCON ||
      _product.get_gen_type() == Product::GenerationType::UNCON_DEG_1 ||
      _product.get_gen_type() == Product::GenerationType::UNCON_SUB) &&
      _product.is_complete()) {
    NodeVector clique;
    for (NodeIt v(_product.get_graph()); v != lemon::INVALID; ++v) {
      clique.push_back(v);
    }
    int size = _product.get_clique_size(clique);
    if (_min_core_size <= size && size <= _max_core_size) {
      _cliques.push_back(clique);
    }
    return false;
  }

  NodeVector order;
  computeDegeneracy(order);

  BitSet mask(_n);

  bool timeout = false;
  for (auto & v : order) {
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

    timeout |= bkPivot(P, D, R, X, S, start, microseconds);
    mask.set(_nodeToBit[v]);

    if (timeout) {
      break;
    }
  }

  return timeout;

}

size_t BronKerbosch::computeDegeneracy() {
  NodeVector order;
  return computeDegeneracy(order);
}

size_t BronKerbosch::computeDegeneracy(NodeVector& order) {
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
    if (!T[i].empty()) {
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

bool BronKerbosch::bkPivot(BitSet P, const BitSet & D,
                           const BitSet & R,
                           BitSet X, const BitSet & S,
                           std::chrono::high_resolution_clock::time_point start,
                           long microseconds) {
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

  std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
  long duration = std::chrono::duration_cast<std::chrono::microseconds>( now - start ).count();

  // Reports maximal c-cliques in P \cup R (but not in X and S)
  BitSet P_cup_X = P | X;
  if (P_cup_X.none()) {
    report(R);
    return false;
  } else {
    if (microseconds > 0 && duration < microseconds) {
      bool timeout = false;
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
          const BitSet &N_v = _bitNeighborhood[v];
          const BitSet &Nc_v = _restrictedBitNeighborhood[v];

          BitSet P_ = P | (D & Nc_v);
          BitSet D_ = D - Nc_v;

          BitSet R_ = R;
          R_[_nodeToBit[v]] = true;

          BitSet X_ = X | (S & Nc_v);
          BitSet S_ = S - Nc_v;

          // report all maximal cliques in ( (P | N[v]) & R) \ (X & N[v]) )
          timeout |= bkPivot(P_ & N_v,
                             D_ & N_v,
                             R_,
                             X_ & N_v,
                             S_ & N_v,
                             start,
                             microseconds);
          P[i] = false;
          X[i] = true;
        }
      }
      return timeout;
    } else {
      return true;
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

