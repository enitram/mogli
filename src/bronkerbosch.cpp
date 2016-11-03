//
// Created by M. Engler on 02/11/16.
//

#include "bronkerbosch.h"

using namespace mogli;

void BronKerbosch::run(SolverType type) {

  BitSet P(_n), R(_n), X(_n);
  P.set();

  bkClassic(P, R, X);

//  switch (type) {
//    case BK_CLASSIC:
//    case BK_PIVOT: {
//      BitSet P(_n), R(_n), X(_n);
//      P.set();
//
//      if (type == BK_CLASSIC)
//        bkClassic(P, R, X);
//      else
//        bkPivot(P, R, X);
//    }
//      break;
//    case BK_PIVOT_DEGENERACY: {
//      NodeVector order;
//      computeDegeneracy(order);
//      bkDegeneracy(order);
//    }
//      break;
//  }

}

void BronKerbosch::bkClassic(BitSet P, BitSet R, BitSet X) {
  assert((P & X).none());
  assert((P & R).none());
  assert((R & X).none());

  // let's print P, R and X
  //std::cout << "P = ";
  //print(P, std::cout);
  //std::cout << ", R = ";
  //print(R, std::cout);
  //std::cout << ", X = ";
  //print(X, std::cout);
  //std::cout << std::endl;

  // Reports maximal cliques in P \cup R (but not in X)
  BitSet P_cup_X = P | X;
  if (P_cup_X.none()) {
    report(R);
  } else {
    for (size_t i = 0; i < P.size(); ++i) {
      if (P[i]) {
        Node v = _bitToNode[i];
        BitSet R_ = R;
        R_[_nodeToBit[v]] = 1;
        // report all maximal cliques in ( (P | N[v]) & R) \ (X & N[v]) )
        bkClassic(P & _bitNeighborhood[v], R_, X & _bitNeighborhood[v]);
        P[i] = 0;
        X[i] = 1;
      }
    }
  }

  /* Invariants:
   * - R is a clique
   * - Each node v in P is adjacent to all nodes in R
   *   (nodes in P are used to extend R, upon usage it's moved to X)
   * - Each node v in X is adjacent to all nodes in R
   *   (maximal cliques containing R \cup v have already been reported)
   */
}

void BronKerbosch::report(const BitSet& R) {
  NodeVector clique;
  for (size_t i = 0; i < R.size(); ++i) {
    if (R[i]) {
//      if (g_verbosity >= VERBOSE_DEBUG)
//      {
//        std::cerr << " " << _g.id(_bitToNode[i]);
//      }
      clique.push_back(_bitToNode[i]);
    }
  }
  _cliques.push_back(clique);
//  if (g_verbosity >= VERBOSE_DEBUG)
//    std::cerr << std::endl;
}

