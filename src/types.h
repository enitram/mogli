//
// Created by M. Engler on 07/12/16.
//

#ifndef MOGLI_TYPES_H
#define MOGLI_TYPES_H

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
  typedef typename Graph::template NodeMap<bool> NodeToBoolMap;

  typedef std::vector<bool> BoolVector;
  typedef std::vector<unsigned short> ShortVector;
  typedef std::vector<int> IntVector;
  typedef std::vector<unsigned long> LongVector;

}

#endif //MOGLI_TYPES_H
