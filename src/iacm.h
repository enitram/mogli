//
// Created by M. Engler on 25/10/16.
//

#ifndef MOGLI_BOOSTERPACK_IACM_H
#define MOGLI_BOOSTERPACK_IACM_H

class IACM {

public:

  std::map<std::string, int> IACM_ELEMENTS;
  std::map<int, std::string> IACM_ELEMENT_NUMBERS;

  std::map<std::string, int> CHEM_ELEMENTS;
  std::map<int, std::string> CHEM_ELEMENT_NUMBERS;

  std::map<std::string, std::string> CPK_COLORS;

  std::map<std::string, int>::const_iterator IACM_ELEMENTS_END;
  std::map<int, std::string>::const_iterator IACM_ELEMENT_NUMBERS_END;

  std::map<std::string, int>::const_iterator CHEM_ELEMENTS_END;
  std::map<int, std::string>::const_iterator CHEM_ELEMENT_NUMBERS_END;

  std::map<std::string, std::string>::const_iterator CPK_COLORS_END;

  std::string UNKNOWN;
  std::string UNKNOWN_COLOR;
  unsigned short UNKNOWN_NUMBER;

private:

  IACM() : IACM_ELEMENTS({
    {"O", 1},
    {"OM", 2},
    {"OA", 3},
    {"OE", 4},
    {"OW", 5},
    {"N", 6},
    {"NT", 7},
    {"NL", 8},
    {"NR", 9},
    {"NZ", 10},
    {"NE", 11},
    {"C", 12},
    {"CH0", 13},
    {"CH1", 14},
    {"CH2", 15},
    {"CH3", 16},
    {"CH4", 17},
    {"CH2r", 18},
    {"CR1", 19},
    {"HC", 20},
    {"H", 21},
    {"DUM", 22},
    {"S", 23},
    {"CU1+", 24},
    {"CU2+", 25},
    {"FE", 26},
    {"ZN2+", 27},
    {"MG2+", 28},
    {"CA2+", 29},
    {"P,SI", 30},
    {"AR", 31},
    {"F", 32},
    {"CL", 33},
    {"BR", 34},
    {"CMet", 35},
    {"OMet", 36},
    {"NA+", 37},
    {"CL-", 38},
    {"CChl", 39},
    {"CLChl", 40},
    {"HChl", 41},
    {"SDmso", 42},
    {"CDmso", 43},
    {"ODmso", 44},
    {"CCl4", 45},
    {"CLCl4", 46},
    {"FTFE", 47},
    {"CTFE", 48},
    {"CHTFE", 49},
    {"OTFE", 50},
    {"CUrea", 51},
    {"OUrea", 52},
    {"NUrea", 53},
    {"CH3p", 54},
    {"I", 55},
    {"CLOpt", 56},
    {"B", 57},
    {"SE", 58},
    {"HS14", 59},
    {"CLAro", 60},
    {"BROpt", 61}
  }), IACM_ELEMENT_NUMBERS({
    {1, "O"},
    {2, "OM"},
    {3, "OA"},
    {4,  "OE"},
    {5, "OW"},
    {6,"N"},
    {7, "NT"},
    {8, "NL"},
    {9, "NR"},
    {10, "NZ"},
    {11, "NE"},
    {12, "C"},
    {13, "CH0"},
    {14, "CH1"},
    {15, "CH2"},
    {16, "CH3"},
    {17, "CH4"},
    {18, "CH2r"},
    {19, "CR1"},
    {20, "HC"},
    {21, "H"},
    {22, "DUM"},
    {23, "S"},
    {24, "CU1+"},
    {25, "CU2+"},
    {26, "FE"},
    {27, "ZN2+"},
    {28, "MG2+"},
    {29, "CA2+"},
    {30, "P,SI"},
    {31, "AR"},
    {32, "F"},
    {33, "CL"},
    {34, "BR"},
    {35, "CMet"},
    {36, "OMet"},
    {37, "NA+"},
    {38, "CL-"},
    {39, "CChl"},
    {40, "CLChl"},
    {41, "HChl"},
    {42, "SDmso"},
    {43, "CDmso"},
    {44, "ODmso"},
    {45, "CCl4"},
    {46, "CLCl4"},
    {47, "FTFE"},
    {48, "CTFE"},
    {49, "CHTFE"},
    {50, "OTFE"},
    {51, "CUrea"},
    {52, "OUrea"},
    {53, "NUrea"},
    {54, "CH3p"},
    {55, "I"},
    {56, "CLOpt"},
    {57, "B"},
    {58, "SE"},
    {59, "HS14"},
    {60, "CLAro"},
    {61, "BROpt"}
  }), CHEM_ELEMENTS({
    {"O", 1},
    {"N", 6},
    {"C", 12},
    {"H", 20},
    {"S", 23},
    {"Cu", 24},
    {"Fe", 26},
    {"Zn", 27},
    {"Mg", 28},
    {"Ca", 29},
    {"P", 30},
    {"Ar", 31},
    {"F", 32},
    {"Cl", 33},
    {"Br", 34},
    {"Na", 37},
    {"B", 57},
    {"Se", 58},
  }), CHEM_ELEMENT_NUMBERS({
    {1, "O"},
    {2, "O"},
    {3, "O"},
    {4, "O"},
    {5, "O"},
    {6, "N"},
    {7, "N"},
    {8, "N"},
    {9, "N"},
    {10, "N"},
    {11, "N"},
    {12, "C"},
    {13, "C"},
    {14, "C"},
    {15, "C"},
    {16, "C"},
    {17, "C"},
    {18, "C"},
    {19, "C"},
    {20, "H"},
    {21, "H"},
    {23, "S"},
    {24, "Cu"},
    {25, "Cu"},
    {26, "Fe"},
    {27, "Zn"},
    {28, "Mg"},
    {29, "Ca"},
    {30, "P"},
    {31, "Ar"},
    {32, "F"},
    {33, "Cl"},
    {34, "Br"},
    {35, "C"},
    {36, "O"},
    {37, "Na"},
    {38, "Cl"},
    {39, "C"},
    {40, "Cl"},
    {41, "H"},
    {42, "S"},
    {43, "C"},
    {44, "O"},
    {45, "C"},
    {46, "Cl"},
    {47, "F"},
    {48, "C"},
    {49, "C"},
    {50, "O"},
    {51, "C"},
    {52, "O"},
    {53, "N"},
    {54, "C"},
    {55, "I"},
    {56, "Cl"},
    {57, "B"},
    {58, "Se"},
    {59, "H"},
    {60, "Cl"},
    {61, "Br"}
  }), CPK_COLORS({
    {"C", "grey54"},
    {"H", "white"},
    {"N", "blue"},
    {"O", "red"},
    {"P", "orange"},
    {"S", "yellow"},
    {"F", "green"},
    {"Cl", "green"},
    {"Br", "red3"},
    {"I", "palevioletred3"},
    {"Ar", "cyan"},
    {"B", "salmon"},
    {"Na", "blueviolet"},
    {"Mg", "green4"},
    {"Fe", "darkorange"}
  }),
     IACM_ELEMENTS_END(IACM_ELEMENTS.end()),
     IACM_ELEMENT_NUMBERS_END(IACM_ELEMENT_NUMBERS.end()),
     CHEM_ELEMENTS_END(CHEM_ELEMENTS.end()),
     CHEM_ELEMENT_NUMBERS_END(CHEM_ELEMENT_NUMBERS.end()),
     CPK_COLORS_END(CPK_COLORS.end()),
     UNKNOWN("*"),
     UNKNOWN_COLOR("pink"),
     UNKNOWN_NUMBER(std::numeric_limits<unsigned short>::max())
  {}

public:

  const std::string get_iacm_element(unsigned short num) {
    auto it = IACM_ELEMENT_NUMBERS.find(num);
    if (it  != IACM_ELEMENT_NUMBERS_END)
      return (*it).second;
    return UNKNOWN;
  }

  const std::string get_chem_element(unsigned short num) {
    auto it = CHEM_ELEMENT_NUMBERS.find(num);
    if (it  != CHEM_ELEMENT_NUMBERS_END)
      return (*it).second;
    return UNKNOWN;
  }

  const std::string get_chem_color(unsigned short num) {
    std::string element = get_chem_element(num);
    auto it = CPK_COLORS.find(element);
    if (it  != CPK_COLORS_END)
      return (*it).second;
    else return UNKNOWN_COLOR;
  }

  const unsigned short get_number(std::string element) {
    auto it = IACM_ELEMENTS.find(element);
    if (it != IACM_ELEMENTS_END)
      return static_cast<unsigned short>((*it).second);
    auto it2 = CHEM_ELEMENTS.find(element);
    if (it2 != CHEM_ELEMENTS_END) {
      return static_cast<unsigned short>((*it2).second);
    }
    return UNKNOWN_NUMBER;
  }

  static IACM &get_default() {
    static IACM instance;
    return instance;
  }

  IACM(IACM const&) = delete;
  void operator=(IACM const&)  = delete;

};


#endif //MOGLI_BOOSTERPACK_IACM_H
