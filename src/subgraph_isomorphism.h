//
// Created by M. Engler on 06/12/16.
//

#ifndef MOGLI_SUBGRAPH_ISOMORPHISM_H
#define MOGLI_SUBGRAPH_ISOMORPHISM_H

#include <sublad.h>
#include <malloc.h>
#include "molecule.h"

namespace mogli {

  Tgraph* translate_graph(const Molecule &mol);
  
  void free_graph(Tgraph* graph);
  
  bool are_subgraph_isomorphic(const Molecule &mol_small, const Molecule &mol_large, int map[]);

  bool are_subgraph_isomorphic(Tgraph* graph_small, Tgraph* graph_large, int map[]);

}

#endif //MOGLI_SUBGRAPH_ISOMORPHISM_H
