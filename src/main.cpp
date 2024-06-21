#include "definitions.hpp"
#include "Nexullance_IT.hpp"
#include "lib_nexullance_IT.hpp"
#include <iostream>
#include <fstream>
#include <string>

#include <pybind11/pybind11.h>
// #include <pybind11/eigen.h>
// #include "mylib.h"
namespace py = pybind11;
// constexpr auto byref = py::return_value_policy::reference_internal;

int test_func() {

    // if (argc < 4) {
    //     std::cerr << "Usage: " << argv[0] << " <input_graph_path>" << " <input_demand_matrix_path>" << "<debug?>"<< std::endl;
    //     return 1;
    // }

    // std::vector<std::pair<int, int>> configs = {{25, 6}};
    std::vector<std::pair<int, int>> configs = {{16, 5}, {25, 6}, {36, 7}, {49, 8}, {64, 9}, {81, 10}, {100, 11}};
    std::vector<std::string> traffics = {"uniform", "shift_1", "shift_half"};

    std::ofstream output_file("IT_boost.csv");
    output_file << "V,D,traffic_pattern,time[s],max_load" << std::endl;
    for (const auto& config : configs) {
        int V = config.first;
        int D = config.second;
        std::string input_graph_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/graph.graphml";
        
        for (const auto& traffic : traffics) {
            std::string input_demand_matrix_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/" + traffic+ ".txt";
            
            double ave_time = 0;
            float ave_max_load = 0;
            // iterate ten times:
            for (int i = 0; i < 10; i++) {
            std::tuple<double, float> result = run_Nexullance_IT_with_paths(input_graph_path, input_demand_matrix_path, 1);
            ave_time += std::get<0>(result);  // Access the first element
            ave_max_load += std::get<1>(result);  // Access the second element
            }
            ave_time/=10;
            ave_max_load/=10;

            // write into the csv file:
            output_file << V << "," << D << "," << traffic << "," << ave_time << "," << ave_max_load << std::endl;
        }
    }
    output_file.close();


    // std::ofstream output_file_2("2+6_IT_boost.csv");
    // output_file_2 << "V,D,traffic_pattern,time[s],max_load" << std::endl;
    // for (const auto& config : configs) {
    //     int V = config.first;
    //     int D = config.second;
    //     std::string input_graph_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/graph.graphml";
        
    //     for (const auto& traffic : traffics) {
    //         std::string input_demand_matrix_path = "/users/ziyzhang/topology-research/nexullance/handover_data/RRG_" + std::to_string(V) + "_" + std::to_string(D) + "/" + traffic+ ".txt";
    //         std::tuple<double, float> result = run_Nexullance_IT(input_graph_path, input_demand_matrix_path, 2);
    //         double time = std::get<0>(result);  // Access the first element
    //         float max_load = std::get<1>(result);  // Access the second element

    //         // write into the csv file:
    //         output_file_2 << V << "," << D << "," << traffic << "," << time << "," << max_load << std::endl;
    //     }
    // }
    // output_file_2.close();
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



// PYBIND11_MODULE(Nexullance_IT, m) { //TODO: move to another cpp file
//     m.doc() = "optional module docstring";

//     m.def("test_func", &test_func, "test_func_cpp_to_py");
//     // py::class_<MyClass>(m, "MyClass")
//     // .def(py::init<double, double, int>())  
//     // .def("run", &MyClass::run, py::call_guard<py::gil_scoped_release>())
//     // .def_readonly("v_data", &MyClass::v_data, byref)
//     // .def_readonly("v_gamma", &MyClass::v_gamma, byref)
//     // ;
// }