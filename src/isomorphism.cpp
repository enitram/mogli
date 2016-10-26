//
// Created by martin on 10/21/16.
//

#include "isomorphism.h"

bool mogli::are_isomorph(Molecule &mol1, Molecule &mol2) {
  if (mol1.get_atom_count() != mol2.get_atom_count())
    return false;
  Canonization c1(mol1);
  Canonization c2(mol2);
  return are_isomorph(c1, c2);
}

bool mogli::are_isomorph(mogli::Canonization &canon1, mogli::Canonization &canon2) {
  const ShortVector& colors1 = canon1.get_colors();
  const ShortVector& colors2 = canon2.get_colors();

  if (colors1.size() != colors2.size())
    return false;

  const LongVector& canonization1 = canon1.get_canonization();
  const LongVector& canonization2 = canon2.get_canonization();

  if (canonization1.size() != canonization2.size())
    return false;

  for (ShortVector::const_iterator i1 = colors1.begin(), i2 = colors2.begin(),
           ie1 = colors1.end(), ie2 = colors2.end();
       i1 != ie1 && i2 != ie2; ++i1, ++i2) {
    if (*i1 != *i2)
      return false;
  }

  for (LongVector::const_iterator i1 = canonization1.begin(), i2 = canonization2.begin(),
           ie1 = canonization1.end(), ie2 = canonization2.end();
       i1 != ie1 && i2 != ie2; ++i1, ++i2) {
    if (*i1 != *i2)
      return false;
  }

  return true;
}
