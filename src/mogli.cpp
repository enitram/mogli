//
// Created by martin on 10/20/16.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "isomorphism.h"

using namespace mogli;

int main() {
  std::ifstream file1("data/min_1.lgf"), file2("data/min_2.lgf"), file3("data/min_3.lgf");
  std::stringstream buffer1, buffer2, buffer3;
  buffer1 << file1.rdbuf();
  buffer2 << file2.rdbuf();
  buffer3 << file3.rdbuf();

  std::string lgf1 = buffer1.str();
  std::string lgf2 = buffer2.str();
  std::string lgf3 = buffer3.str();

  Molecule mol1, mol2, mol3;

  mol1.readLGF(lgf1);
  mol2.readLGF(lgf2);
  mol3.readLGF(lgf3);

  std::cout << areIsomorph(mol1, mol2) << std::endl;
  std::cout << areIsomorph(mol1, mol1) << std::endl;
  std::cout << areIsomorph(mol2, mol2) << std::endl;
  std::cout << areIsomorph(mol2, mol3) << std::endl;
  std::cout << areIsomorph(mol3, mol2) << std::endl;

  return 0;
}

