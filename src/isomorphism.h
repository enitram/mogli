//
// Created by M. Engler on 10/21/16.
//

#ifndef MOGLI_ISOMORPHISM_H
#define MOGLI_ISOMORPHISM_H

#include "canonization.h"

namespace mogli {

  bool are_isomorphic(Molecule &mol1, Molecule &mol2);

  bool are_isomorphic(Canonization &canon1, Canonization &canon2);

}

#endif //MISO_ISOMORPHISM_H
