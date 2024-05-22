#include <iostream>
#include <fstream> 
#include <boost/graph/graphml.hpp>
#include <boost/property_map/dynamic_property_map.hpp>  // For dynamic_properties
#include <boost/property_map/property_map.hpp>          // For property

#include "definitions.hpp"

Graph read_graphml(std::string graphmlFile, bool print) {

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
    // Read the GraphML file into the graph
    read_graphml(inFile, g, dp);
    inFile.close();

    // Assign all edge weights as 1.0
    auto edge_pair = edges(g);
    for (auto it = edge_pair.first; it != edge_pair.second; ++it) {
        put(boost::edge_weight, g, *it, 1.0f);
    }

    if (print){
        std::cout << "Graph contents:" << std::endl;
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

        // Iterate through the edges and print their source, destination, and weight
        for (auto it = edge_pair.first; it != edge_pair.second; ++it) {
            auto src = source(*it, g);
            auto dest = target(*it, g);
            auto weight = get(boost::edge_weight, g, *it);
            std::cout << "Edge from vertex " << src << " to vertex " << dest
                    << " with weight " << weight << std::endl;
    }

    }
    return g;
}


void read_matrix(std::string filename, bool print, float** result_matrix, size_t dim){
    // read the matrix from the file to result_matrix
    std::ifstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (file.good()) {
                    std::string value_str = "";
                    if (j < dim - 1) {
                        std::getline(file, value_str, ',');
                    } else {
                        std::getline(file, value_str, '\n');
                    }
                    // assert that the matrix has the correct dimensions
                    assert(value_str.size() > 0 && "is the matrix dimension correct?");
                    result_matrix[i][j] = std::stof(value_str);
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Failed to open file: " << filename << std::endl;
    }

    if (print) {
        std::cout << "Matrix contents:" << std::endl;
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                std::cout << result_matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
}

