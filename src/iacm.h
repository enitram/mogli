//
// Created by martin on 25/10/16.
//

#ifndef MOGLI_BOOSTERPACK_IACM_H
#define MOGLI_BOOSTERPACK_IACM_H

class IACM {

public:

  std::map<std::string, unsigned short> IACM_ELEMENTS;
  std::map<unsigned short, std::string> IACM_ELEMENT_NUMBERS;

  std::map<std::string, unsigned short> CHEM_ELEMENTS;
  std::map<unsigned short, std::string> CHEM_ELEMENT_NUMBERS;

  std::map<std::string, unsigned short>::const_iterator IACM_ELEMENTS_END;
  std::map<unsigned short, std::string>::const_iterator IACM_ELEMENT_NUMBERS_END;

  std::map<std::string, unsigned short>::const_iterator CHEM_ELEMENTS_END;
  std::map<unsigned short, std::string>::const_iterator CHEM_ELEMENT_NUMBERS_END;

  std::string UNKNOWN;
  unsigned short UNKNOWN_NUMBER;

private:

  IACM() : IACM_ELEMENTS({
     {"O", static_cast<unsigned short>(1)},
     {"OM", static_cast<unsigned short>(2)},
     {"OA", static_cast<unsigned short>(3)},
     {"OE", static_cast<unsigned short>(4)},
     {"OW", static_cast<unsigned short>(5)},
     {"N", static_cast<unsigned short>(6)},
     {"NT", static_cast<unsigned short>(7)},
     {"NL", static_cast<unsigned short>(8)},
     {"NR", static_cast<unsigned short>(9)},
     {"NZ", static_cast<unsigned short>(10)},
     {"NE", static_cast<unsigned short>(11)},
     {"C", static_cast<unsigned short>(12)},
     {"CH0", static_cast<unsigned short>(13)},
     {"CH1", static_cast<unsigned short>(14)},
     {"CH2", static_cast<unsigned short>(15)},
     {"CH3", static_cast<unsigned short>(16)},
     {"CH4", static_cast<unsigned short>(17)},
     {"CH2r", static_cast<unsigned short>(18)},
     {"CR1", static_cast<unsigned short>(19)},
     {"HC", static_cast<unsigned short>(20)},
     {"H", static_cast<unsigned short>(21)},
     {"DUM", static_cast<unsigned short>(22)},
     {"S", static_cast<unsigned short>(23)},
     {"CU1+", static_cast<unsigned short>(24)},
     {"CU2+", static_cast<unsigned short>(25)},
     {"FE", static_cast<unsigned short>(26)},
     {"ZN2+", static_cast<unsigned short>(27)},
     {"MG2+", static_cast<unsigned short>(28)},
     {"CA2+", static_cast<unsigned short>(29)},
     {"P,SI", static_cast<unsigned short>(30)},
     {"AR", static_cast<unsigned short>(31)},
     {"F", static_cast<unsigned short>(32)},
     {"CL", static_cast<unsigned short>(33)},
     {"BR", static_cast<unsigned short>(34)},
     {"CMet", static_cast<unsigned short>(35)},
     {"OMet", static_cast<unsigned short>(36)},
     {"NA+", static_cast<unsigned short>(37)},
     {"CL-", static_cast<unsigned short>(38)},
     {"CChl", static_cast<unsigned short>(39)},
     {"CLChl", static_cast<unsigned short>(40)},
     {"HChl", static_cast<unsigned short>(41)},
     {"SDmso", static_cast<unsigned short>(42)},
     {"CDmso", static_cast<unsigned short>(43)},
     {"ODmso", static_cast<unsigned short>(44)},
     {"CCl4", static_cast<unsigned short>(45)},
     {"CLCl4", static_cast<unsigned short>(46)},
     {"FTFE", static_cast<unsigned short>(47)},
     {"CTFE", static_cast<unsigned short>(48)},
     {"CHTFE", static_cast<unsigned short>(49)},
     {"OTFE", static_cast<unsigned short>(50)},
     {"CUrea", static_cast<unsigned short>(51)},
     {"OUrea", static_cast<unsigned short>(52)},
     {"NUrea", static_cast<unsigned short>(53)},
     {"CH3p", static_cast<unsigned short>(54)},
     {"I", static_cast<unsigned short>(55)},
     {"CLOpt", static_cast<unsigned short>(56)},
     {"B", static_cast<unsigned short>(57)},
     {"SE", static_cast<unsigned short>(58)},
     {"HS14", static_cast<unsigned short>(59)},
     {"CLAro", static_cast<unsigned short>(60)},
     {"BROpt", static_cast<unsigned short>(61)}
  }), IACM_ELEMENT_NUMBERS({
    {static_cast<unsigned short>(1), "O"},
    {static_cast<unsigned short>(2), "OM"},
    {static_cast<unsigned short>(3), "OA"},
    {static_cast<unsigned short>(4),  "OE"},
    {static_cast<unsigned short>(5), "OW"},
    {static_cast<unsigned short>(6),"N"},
    {static_cast<unsigned short>(7), "NT"},
    {static_cast<unsigned short>(8), "NL"},
    {static_cast<unsigned short>(9), "NR"},
    {static_cast<unsigned short>(10), "NZ"},
    {static_cast<unsigned short>(11), "NE"},
    {static_cast<unsigned short>(12), "C"},
    {static_cast<unsigned short>(13), "CH0"},
    {static_cast<unsigned short>(14), "CH1"},
    {static_cast<unsigned short>(15), "CH2"},
    {static_cast<unsigned short>(16), "CH3"},
    {static_cast<unsigned short>(17), "CH4"},
    {static_cast<unsigned short>(18), "CH2r"},
    {static_cast<unsigned short>(19), "CR1"},
    {static_cast<unsigned short>(20), "HC"},
    {static_cast<unsigned short>(21), "H"},
    {static_cast<unsigned short>(22), "DUM"},
    {static_cast<unsigned short>(23), "S"},
    {static_cast<unsigned short>(24), "CU1+"},
    {static_cast<unsigned short>(25), "CU2+"},
    {static_cast<unsigned short>(26), "FE"},
    {static_cast<unsigned short>(27), "ZN2+"},
    {static_cast<unsigned short>(28), "MG2+"},
    {static_cast<unsigned short>(29), "CA2+"},
    {static_cast<unsigned short>(30), "P,SI"},
    {static_cast<unsigned short>(31), "AR"},
    {static_cast<unsigned short>(32), "F"},
    {static_cast<unsigned short>(33), "CL"},
    {static_cast<unsigned short>(34), "BR"},
    {static_cast<unsigned short>(35), "CMet"},
    {static_cast<unsigned short>(36), "OMet"},
    {static_cast<unsigned short>(37), "NA+"},
    {static_cast<unsigned short>(38), "CL-"},
    {static_cast<unsigned short>(39), "CChl"},
    {static_cast<unsigned short>(40), "CLChl"},
    {static_cast<unsigned short>(41), "HChl"},
    {static_cast<unsigned short>(42), "SDmso"},
    {static_cast<unsigned short>(43), "CDmso"},
    {static_cast<unsigned short>(44), "ODmso"},
    {static_cast<unsigned short>(45), "CCl4"},
    {static_cast<unsigned short>(46), "CLCl4"},
    {static_cast<unsigned short>(47), "FTFE"},
    {static_cast<unsigned short>(48), "CTFE"},
    {static_cast<unsigned short>(49), "CHTFE"},
    {static_cast<unsigned short>(50), "OTFE"},
    {static_cast<unsigned short>(51), "CUrea"},
    {static_cast<unsigned short>(52), "OUrea"},
    {static_cast<unsigned short>(53), "NUrea"},
    {static_cast<unsigned short>(54), "CH3p"},
    {static_cast<unsigned short>(55), "I"},
    {static_cast<unsigned short>(56), "CLOpt"},
    {static_cast<unsigned short>(57), "B"},
    {static_cast<unsigned short>(58), "SE"},
    {static_cast<unsigned short>(59), "HS14"},
    {static_cast<unsigned short>(60), "CLAro"},
    {static_cast<unsigned short>(61), "BROpt"}
  }), CHEM_ELEMENTS({
    {"O", static_cast<unsigned short>(1)},
    {"N", static_cast<unsigned short>(6)},
    {"C", static_cast<unsigned short>(12)},
    {"H", static_cast<unsigned short>(20)},
    {"S", static_cast<unsigned short>(23)}
  }), CHEM_ELEMENT_NUMBERS({
     {static_cast<unsigned short>(1), "O"},
     {static_cast<unsigned short>(2), "O"},
     {static_cast<unsigned short>(3), "O"},
     {static_cast<unsigned short>(4), "O"},
     {static_cast<unsigned short>(5), "O"},
     {static_cast<unsigned short>(6), "N"},
     {static_cast<unsigned short>(7), "N"},
     {static_cast<unsigned short>(8), "N"},
     {static_cast<unsigned short>(9), "N"},
     {static_cast<unsigned short>(10), "N"},
     {static_cast<unsigned short>(11), "N"},
     {static_cast<unsigned short>(12), "C"},
     {static_cast<unsigned short>(13), "C"},
     {static_cast<unsigned short>(14), "C"},
     {static_cast<unsigned short>(15), "C"},
     {static_cast<unsigned short>(16), "C"},
     {static_cast<unsigned short>(17), "C"},
     {static_cast<unsigned short>(18), "C"},
     {static_cast<unsigned short>(19), "C"},
     {static_cast<unsigned short>(20), "H"},
     {static_cast<unsigned short>(21), "H"},
     {static_cast<unsigned short>(23), "S"},
     {static_cast<unsigned short>(42), "S"}
   }),
     IACM_ELEMENTS_END(IACM_ELEMENTS.end()),
     IACM_ELEMENT_NUMBERS_END(IACM_ELEMENT_NUMBERS.end()),
     CHEM_ELEMENTS_END(CHEM_ELEMENTS.end()),
     CHEM_ELEMENT_NUMBERS_END(CHEM_ELEMENT_NUMBERS.end()),
     UNKNOWN("?"),
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

  const unsigned short get_number(std::string element) {
    auto it = IACM_ELEMENTS.find(element);
    if (it != IACM_ELEMENTS_END)
      return (*it).second;
    auto it2 = CHEM_ELEMENTS.find(element);
    if (it2 != CHEM_ELEMENTS_END) {
      return (*it2).second;
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
