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
  typedef std::vector<std::string> StringVector;
  typedef typename Graph::template NodeMap<bool> NodeToBoolMap;

}

#endif //MOGLI_TYPES_H
