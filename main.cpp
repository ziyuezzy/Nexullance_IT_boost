// #include <iostream>
// #include <fstream>  // For file stream operations
// #include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/graphml.hpp>
// #include <boost/property_map/dynamic_property_map.hpp>  // For dynamic_properties
// #include <boost/property_map/property_map.hpp>          // For property

// using namespace boost;

#include "definitions.hpp"

// Define a directed graph with integer vertex and edge properties
typedef adjacency_list<vecS, vecS, directedS,
                       property<vertex_index_t, int>,
                       property<edge_index_t, int>
                      > Graph;

int main() {
    read_and_print_graphml("/users/ziyzhang/topology-research/graphml_data/RRG_18_5.graphml", true);
    return 0;
}
