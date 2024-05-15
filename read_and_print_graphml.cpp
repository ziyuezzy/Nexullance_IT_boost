#include <iostream>
#include <fstream>  // For file stream operations
// #include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>  // For dynamic_properties
#include <boost/property_map/property_map.hpp>          // For property

// using namespace boost;

// // Define a directed graph with integer vertex and edge properties
// typedef adjacency_list<vecS, vecS, directedS,
//                        property<vertex_index_t, int>,
//                        property<edge_index_t, int>
//                       > Graph;

#include "definitions.hpp"

Graph read_and_print_graphml(std::string graphmlFile, bool print) {

    // Open the GraphML file for reading
    std::ifstream inFile(graphmlFile);
    if (!inFile.is_open()) {
        std::cerr << "Failed to open file: " << graphmlFile << std::endl;
        return 1;
    }
    // Create a graph object
    Graph g;
    // Setup dynamic properties to read vertex and edge properties
    dynamic_properties dp;
    dp.property("vertex_index", get(vertex_index, g));  // Read vertex indices
    dp.property("edge_index", get(edge_index, g));      // Read edge indices
    // Read the GraphML file into the graph
    read_graphml(inFile, g, dp);

    if (print){
        // Get the range of vertices
        auto vertexRange = boost::vertices(g);

        // Iterate over all vertices and print them
        std::cout << "Vertices:" << std::endl;
        for (auto it = vertexRange.first; it != vertexRange.second; ++it) {
            std::cout << *it << std::endl;
        }

        // Get the range of edges
        auto edgeRange = boost::edges(g);

        // Print out all edges in the format (int, int)
        std::cout << "Edges:" << std::endl;
        for (auto it = edgeRange.first; it != edgeRange.second; ++it) {
            auto source = boost::source(*it, g);
            auto target = boost::target(*it, g);
            std::cout << "(" << source << ", " << target << ")" << std::endl;
        }
    }
    return g;
}
