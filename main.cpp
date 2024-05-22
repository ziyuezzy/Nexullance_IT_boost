#include "definitions.hpp"
#include "Nexullance_IT.hpp"
#include <iostream>
#include <fstream>
#include <string>


int main(int argc, char* argv[]) {

    // if (argc < 4) {
    //     std::cerr << "Usage: " << argv[0] << " <input_graph_path>" << " <input_demand_matrix_path>" << "<debug?>"<< std::endl;
    //     return 1;
    // }

    // std::vector<std::pair<int, int>> configs = {{25, 6}};
    std::vector<std::pair<int, int>> configs = {{16, 5}, {25, 6}, {36, 7}, {49, 8}, {64, 9}, {81, 10}, {100, 11}};
    std::vector<std::string> traffics = {"uniform", "shift_1", "shift_half"};

    std::ofstream output_file("1+6_IT_boost.csv");
    output_file << "V,D,traffic_pattern,time[s],max_load" << std::endl;
    for (const auto& config : configs) {
        int V = config.first;
        int D = config.second;
        std::string input_graph_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/graph.graphml";
        
        for (const auto& traffic : traffics) {
            std::string input_demand_matrix_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/" + traffic+ ".txt";
            std::tuple<double, float> result = run_Nexullance_IT(input_graph_path, input_demand_matrix_path, 1);
            double time = std::get<0>(result);  // Access the first element
            float max_load = std::get<1>(result);  // Access the second element

            // write into the csv file:
            output_file << V << "," << D << "," << traffic << "," << time << "," << max_load << std::endl;
        }
    }
    output_file.close();


    std::ofstream output_file_2("2+6_IT_boost.csv");
    output_file_2 << "V,D,traffic_pattern,time[s],max_load" << std::endl;
    for (const auto& config : configs) {
        int V = config.first;
        int D = config.second;
        std::string input_graph_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/graph.graphml";
        
        for (const auto& traffic : traffics) {
            std::string input_demand_matrix_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/" + traffic+ ".txt";
            std::tuple<double, float> result = run_Nexullance_IT(input_graph_path, input_demand_matrix_path, 2);
            double time = std::get<0>(result);  // Access the first element
            float max_load = std::get<1>(result);  // Access the second element

            // write into the csv file:
            output_file_2 << V << "," << D << "," << traffic << "," << time << "," << max_load << std::endl;
        }
    }
    output_file_2.close();
    return 0;





    // // testing shortest paths implementation (~1000x faster the networkx with -O3 optimization)
    // if (argc < 2) {
    //     std::cerr << "Usage: " << argv[0] << " <input_graph_path>" << std::endl;
    //     return 1;
    // }
    // std::string input_graph_path = argv[1];
    // test_shortest_paths(input_graph_path);

    // return 0;
}