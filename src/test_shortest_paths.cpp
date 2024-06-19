#include "definitions.hpp"
// #include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <string>
#include <chrono>

void test_shortest_paths(std::string input_graph_path) {

    // Graph G = read_graphml("/users/ziyzhang/topology-research/nexullance/handover_data/test_ring/graph.graphml", true);
    Graph G = read_graphml(input_graph_path, false);

    property_map< Graph, edge_weight_t >::type weightmap = get(edge_weight, G);

    bool print_results = true;

    // Define the source vertex (e.g., vertex 0)
    Vertex source = 10;

    std::vector<std::vector<Vertex>>** all_paths_all_s_d = new std::vector<std::vector<Vertex>>*[num_vertices(G)];
    for (int i = 0; i < num_vertices(G); ++i) {
        all_paths_all_s_d[i] = new std::vector<std::vector<Vertex>>[num_vertices(G)];
    }

    auto start = std::chrono::high_resolution_clock::now();
    compute_all_shortest_paths_all_s_d(G, all_paths_all_s_d, weightmap);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Runtime for all s-d dijkstra: " << elapsed.count() << " seconds" << std::endl;

    // // Print all paths

    std::vector<std::vector<Vertex>> all_paths_from_s[num_vertices(G)];

    start = std::chrono::high_resolution_clock::now();
    compute_all_shortest_paths_single_source(G, source, all_paths_from_s, weightmap);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Runtime for single s dijkstra: " << elapsed.count() << " seconds" << std::endl;

    // Print all paths
    if (print_results){
        std::cout << "results all shortest paths:" << std::endl;
        for (size_t i = 0; i < num_vertices(G); ++i) {
            std::cout << "All shortest paths from vertex " << source << " to vertex " << i << " are:" << std::endl;
            for (auto path : all_paths_from_s[i]) {
                std::cout << "Path: ";
                for (auto v : path) {
                    std::cout << v << " ";
                }
                std::cout << std::endl;
            }
        }
    }

    Vertex destination = 8;
    std::vector<std::vector<Vertex>> all_paths;

    start = std::chrono::high_resolution_clock::now();
    compute_all_shortest_paths_single_s_d(G, source, destination, all_paths, weightmap);
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Runtime for single s-d dijkstra: " << elapsed.count() << " seconds" << std::endl;

    if (print_results){
        std::cout << "All shortest paths from vertex " << source << " to vertex " << destination << " are:" << std::endl;
        for (auto path : all_paths) {
            std::cout << "Path: ";
            for (auto v : path) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }
    }

    return;
}
