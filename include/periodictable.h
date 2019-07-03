////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    mogli - molecular graph library                                                                                 //
//                                                                                                                    //
//    Copyright (C) 2016-2019  Martin S. Engler                                                                       //
//                                                                                                                    //
//    This program is free software: you can redistribute it and/or modify                                            //
//    it under the terms of the GNU Lesser General Public License as published                                        //
//    by the Free Software Foundation, either version 3 of the License, or                                            //
//    (at your option) any later version.                                                                             //
//                                                                                                                    //
//    This program is distributed in the hope that it will be useful,                                                 //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                  //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                                    //
//    GNU General Public License for more details.                                                                    //
//                                                                                                                    //
//    You should have received a copy of the GNU Lesser General Public License                                        //
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.                                          //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MOGLI_PERIODICTABLE_H
#define MOGLI_PERIODICTABLE_H

namespace mogli {

  class PeriodicTable {


  private:

    std::map<std::string, unsigned short> ELEMENTS;
    std::map<unsigned short, std::string> ELEMENT_NUMBERS;
    std::map<unsigned short, std::string> COLORS;
    std::string UNKNOWN;
    std::string UNKNOWN_COLOR;
    unsigned short UNKNOWN_NUMBER;

  public:

    PeriodicTable() :
        ELEMENTS(),
        ELEMENT_NUMBERS(),
        COLORS(),
        UNKNOWN("*"),
        UNKNOWN_COLOR("pink"),
        UNKNOWN_NUMBER(std::numeric_limits<unsigned short>::max()) {}

    PeriodicTable &add_uncolored(unsigned short num, std::string name) {
      ELEMENTS[name] = num;
      ELEMENT_NUMBERS[num] = name;
      return *this;
    }

    PeriodicTable &add(unsigned short num, std::string name, std::string color) {
      ELEMENTS[name] = num;
      ELEMENT_NUMBERS[num] = name;
      COLORS[num] = color;
      return *this;
    }

    const std::string get_element(unsigned short num) {
      auto it = ELEMENT_NUMBERS.find(num);
      if (it != ELEMENT_NUMBERS.end())
        return (*it).second;
      return UNKNOWN;
    }

    const std::string get_color(unsigned short num) {
      auto it = COLORS.find(num);
      if (it != COLORS.end())
        return (*it).second;
      else return UNKNOWN_COLOR;
    }

    const unsigned short get_number(std::string element) {
      auto it = ELEMENTS.find(element);
      if (it != ELEMENTS.end())
        return static_cast<unsigned short>((*it).second);
      else return UNKNOWN_NUMBER;
    }

    static PeriodicTable &get_default() {
      static PeriodicTable _instance;
      if (_instance.ELEMENTS.empty()) {
        //    O
        _instance.add(1, "O", "red").add(2, "OM", "red").add(3, "OA", "red").add(4, "OE", "red").add(5, "OW", "red")
            .add(36, "OMet", "red").add(44, "ODmso", "red").add(50, "OTFE", "red").add(52, "OUrea", "red")
            .add(62, "OEOpt", "red").add(68, "OAlc", "red")
                //    N
            .add(6, "N", "blue").add(7, "NT", "blue").add(8, "NL", "blue").add(9, "NR", "blue").add(10, "NZ", "blue")
            .add(11, "NE", "blue").add(53, "NUrea", "blue").add(63, "NOpt", "blue").add(66, "NPri", "blue")
            .add(67, "NTer", "blue")
                //    C
            .add(12, "C", "grey54").add(13, "CH0", "grey54").add(14, "CH1", "grey54").add(15, "CH2", "grey54")
            .add(16, "CH3", "grey54").add(17, "CH4", "grey54").add(18, "CH2r", "grey54").add(19, "CR1", "grey54")
            .add(35, "CMet", "grey54").add(39, "CChl", "grey54").add(43, "CDmso", "grey54").add(45, "CCl4", "grey54")
            .add(48, "CTFE", "grey54").add(49, "CHTFE", "grey54").add(51, "CUrea", "grey54").add(54, "CH3p", "grey54")
            .add(64, "CAro", "grey54").add(65, "CPos", "grey54")
                //    H
            .add(20, "HC", "white").add(21, "H", "white").add(41, "HChl", "white").add(59, "HS14", "white")
                //    S
            .add(23, "S", "yellow").add(42, "SDmso", "yellow")
                //    Cu
            .add(24, "CU1+", "salmon").add(25, "CU2+", "salmon")
                //    Fe
            .add(26, "FE", "darkorange")
                //    Zn
            .add(27, "ZN2+", "salmon")
                //    Mg
            .add(28, "MG2+", "green4")
                //    Ca
            .add(29, "CA2+", "green4")
                //    P
            .add(30, "P,SI", "orange").add(69, "P", "orange")
                //    Ar
            .add(31, "AR", "lightblue")
                //    F
            .add(32, "F", "green").add(47, "FTFE", "green")
                //    Cl
            .add(33, "CL", "green").add(38, "CL-", "green").add(40, "CLChl", "green").add(46, "CLCl4", "green")
            .add(56, "CLOpt", "green").add(60, "CLAro", "green")
                //    Br
            .add(34, "BR", "red3").add(61, "BROpt", "red3")
                //    Na
            .add(37, "NA+", "blueviolet")
                //    I
            .add(55, "I", "palevioletred3")
                //    B
            .add(57, "B", "salmon")
                //    Se
            .add(58, "SE", "goldenrod1")
                // Si
            .add(70, "Si", "orange")
                //    DUMMY
            .add_uncolored(22, "DUM");
      }
      return _instance;
    }

  };

}


#endif //MOGLI_PERIODICTABLE_H
