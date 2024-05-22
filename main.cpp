#include "definitions.hpp"
#include "Nexullance_IT.hpp"
#include <iostream>
#include <string>

#include <chrono>

int main(int argc, char* argv[]) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_graph_path>" << " <input_demand_matrix_path>" << "<debug?>"<< std::endl;
        return 1;
    }

    bool debug = (std::string(argv[3]) == "debug") ? true : false;


    Graph G = read_graphml(argv[1], debug);

    int num_routers=boost::num_vertices(G);
    // reading matrix from file
    float** matrix = new float*[num_routers];
    for (int i = 0; i < num_routers; i++) {
        matrix[i] = new float[num_routers];
    }
    read_matrix(argv[2], debug, matrix, num_routers);

    const float Cap_link = 10;

    NexullanceIT nexu_it = NexullanceIT(G, const_cast<const float**>(matrix), Cap_link, debug);
    
    auto start = std::chrono::high_resolution_clock::now();
    nexu_it.optimize(1, 1.0, 1.0, 6, 3.0, 7.0, 16);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Runtime for Nexullance_IT: " << elapsed.count() << " seconds" << std::endl;

    std::cout << "Nexullance_IT solution found, max link load = "<< nexu_it.result_max_loads_step_2.back() << std::endl;

    // Deallocate the memory for matrix
    for (int i = 0; i < 16; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;


    // // testing shortest paths implementation (~1000x faster the networkx with -O3 optimization)
    // if (argc < 2) {
    //     std::cerr << "Usage: " << argv[0] << " <input_graph_path>" << std::endl;
    //     return 1;
    // }
    // std::string input_graph_path = argv[1];
    // test_shortest_paths(input_graph_path);

    // return 0;
}