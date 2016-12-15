//
// Created by M. Engler on 07/12/16.
//

#include "sublad.h"

c_bool compatibleVertexLabels(int l1, int l2){
  // return true iff a vertex labelled with l1 may be matched with a vertex labelled with l2
  return (l1 == l2);
}

c_bool compatibleEdgeLabels(int l1, int l2){
  // return true iff an edge labelled with l1 may be matched with an edge labelled with l2
  return (l1 == l2);
}
