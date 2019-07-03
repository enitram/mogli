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

#ifndef MOGLI_TYPES_H
#define MOGLI_TYPES_H

#include <lemon/core.h>
#include <lemon/list_graph.h>
#include <map>
#include <vector>
#include <variant>

namespace mogli {

  typedef lemon::ListGraph Graph;
  typedef Graph::Node Node;
  typedef Graph::NodeIt NodeIt;
  typedef Graph::Edge Edge;
  typedef Graph::EdgeIt EdgeIt;
  typedef Graph::IncEdgeIt IncEdgeIt;
  typedef std::vector<Node> NodeVector;

  typedef std::set<int> IntSet;
  typedef std::map<int, int> IntToIntMap;
  typedef std::vector<IntToIntMap> IntToIntMapVector;
  typedef std::vector<std::string> StringVector;
  typedef std::set<std::string> StringSet;
  typedef typename Graph::template NodeMap<bool> NodeToBoolMap;

  typedef std::vector<bool> BoolVector;
  typedef std::vector<unsigned short> ShortVector;
  typedef std::vector<int> IntVector;
  typedef std::vector<unsigned long> LongVector;

  typedef std::variant<bool, int, double, std::string> Any;

  template <typename K, typename V>
  struct UniquePtrMap {
    typedef std::map<K, std::unique_ptr<V>> type;
  };

  template <typename T>
  struct SharedPtrVector {
    typedef std::vector<std::shared_ptr<T>> type;
  };

}

#endif //MOGLI_TYPES_H
