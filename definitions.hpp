#ifndef GRAPH_DEFINITIONS_HPP
#define GRAPH_DEFINITIONS_HPP

#include <boost/graph/adjacency_list.hpp>
using namespace boost;

// Define a directed graph with integer vertex and edge properties
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                               boost::property<boost::vertex_index_t, int>,
                               boost::property<boost::edge_index_t, int>
                              > Graph;

Graph read_and_print_graphml(std::string graphmlFile, bool print);

#endif // GRAPH_DEFINITIONS_HPP
