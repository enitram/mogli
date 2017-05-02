
#include "molecule.h"
#include "canonization.h"

using namespace mogli;

const bool Molecule::is_isomorphic(Molecule &other) const {
  if (_atom_count != other.get_atom_count())
        return false;
      Canonization c1(*this);
      Canonization c2(other);
      return c1.is_isomorphic(c2);
}