//
// Created by M. Engler on 06/12/16.
//

#ifndef MOGLI_SUBGRAPH_ISOMORPHISM_H
#define MOGLI_SUBGRAPH_ISOMORPHISM_H

#include <sublad.h>
#include <malloc.h>
#include "molecule.h"

namespace mogli {

  Tgraph* translate_graph(const Molecule &mol, IntVector &node_ids);
  
  void free_graph(Tgraph* graph);

  void translate_maps(const IntVector &node_ids_small, const IntVector &node_ids_large,
                      const int in_iso_map[], IntToIntMap &out_iso_map);

  bool are_subgraph_isomorphic(const Molecule &mol_small, const Molecule &mol_large,
                               IntToIntMap isomorphism_map);

  bool are_subgraph_isomorphic(Tgraph* graph_small, Tgraph* graph_large, int isomorphism_map[]);

}

#endif //MOGLI_SUBGRAPH_ISOMORPHISM_H
