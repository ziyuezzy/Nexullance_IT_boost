#include "lib_nexullance_IT.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

std::tuple<double, float> run_Nexullance_IT_with_paths(std::string input_graph_path, std::string input_matrix_path, int num_step_1){
    bool debug = false;

    const float Cap_link = 10;
    Graph G = read_graphml(input_graph_path, debug);
    int num_routers=boost::num_vertices(G);
    // reading matrix from file
    float** matrix = new float*[num_routers];
    for (int i = 0; i < num_routers; i++) {
        matrix[i] = new float[num_routers];
    }
    read_matrix(input_matrix_path, debug, matrix, num_routers);


    Nexullance_IT nexu_it = Nexullance_IT(G, const_cast<const float**>(matrix), Cap_link, debug);
    
    auto start = std::chrono::high_resolution_clock::now();
    nexu_it.optimize(num_step_1, 1.0, 1.0, 6, 0.1, 7.0, num_routers);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // // Deallocate the memory for matrix 
    // for (int i = 0; i < num_routers; i++) {
    //     delete[] matrix[i];
    // }
    // delete[] matrix;

    return std::make_tuple(elapsed.count(), nexu_it.result_max_loads_step_2.back());
}

std::tuple<double, float> run_Nexullance_IT(int V, Eigen::MatrixX2i arcs, Eigen::MatrixXf input_matrix, int num_step_1=1){
    bool debug = false;

    const float Cap_link = 10;
    Graph G = read_graph_from_arcs(V, arcs, debug);
    int num_routers=boost::num_vertices(G);


    // converting the Eigen matrix to a float** matrix
    float** matrix = new float*[num_routers];
    for (int i = 0; i < num_routers; i++) {
        matrix[i] = new float[num_routers];
    }
    assert(input_matrix.rows() == num_routers && input_matrix.cols() == num_routers);
    for (int i = 0; i < num_routers; ++i) {
        for (int j = 0; j < num_routers; ++j) {
            matrix[i][j] = input_matrix(i, j);
        }
    }
    //===========

    Nexullance_IT nexu_it = Nexullance_IT(G, const_cast<const float**>(matrix), Cap_link, debug);
    
    auto start = std::chrono::high_resolution_clock::now();
    nexu_it.optimize(num_step_1, 1.0, 1.0, 6, 0.1, 7.0, num_routers);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // // Deallocate the memory for matrix
    // for (int i = 0; i < num_routers; i++) {
    //     delete[] matrix[i];
    // }
    // delete[] matrix;

    return std::make_tuple(elapsed.count(), nexu_it.result_max_loads_step_2.back());
}

Nexullance_IT_interface::Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, Eigen::MatrixXf input_matrix, bool debug): _V(V) {
    const float Cap_link = 10;
    Graph G = read_graph_from_arcs(V, arcs, false);
    int num_routers=boost::num_vertices(G);
    assert(V == num_routers);
    // converting the Eigen matrix to a float** matrix
    float** matrix = new float*[num_routers];
    for (int i = 0; i < num_routers; i++) {
        matrix[i] = new float[num_routers];
    }
    assert(input_matrix.rows() == num_routers && input_matrix.cols() == num_routers);
    for (int i = 0; i < num_routers; ++i) {
        for (int j = 0; j < num_routers; ++j) {
            matrix[i][j] = input_matrix(i, j);
        }
    }
    //===========
    nexu_it = new Nexullance_IT(G, const_cast<const float**>(matrix), Cap_link, debug);
}

Nexullance_IT_interface::~Nexullance_IT_interface(){
    delete nexu_it;
}

void Nexullance_IT_interface::set_parameters(float _alpha, float _beta){
    alpha = _alpha;
    beta = _beta;
}

double Nexullance_IT_interface::run(){
    auto start = std::chrono::high_resolution_clock::now();
    nexu_it->optimize(1, 1.0, 1.0, 6, alpha, beta, _V);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    return elapsed.count();
}

float Nexullance_IT_interface::get_max_link_load(){
    return nexu_it->result_max_loads_step_2.back(); 
}

Nexullance_IT::result_routing_table Nexullance_IT_interface::get_routing_table(){
    return nexu_it->get_routing_table();
}

size_t Nexullance_IT_interface::get_num_attempts_step_2(){
    return nexu_it->num_attempts_step_2;
}


namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(Nexullance_IT_cpp, m) {
    m.doc() = "calling Nexullance IT";

    m.def("run_Nexullance_IT", &run_Nexullance_IT, "test_func_directly_calling_nexullance_IT");
    m.def("run_Nexullance_IT_with_paths", &run_Nexullance_IT_with_paths, "test_func_directly_calling_nexullance_IT");
    // py::class_<MyClass>(m, "MyClass")
    // .def(py::init<double, double, int>())  
    // .def("run", &MyClass::run, py::call_guard<py::gil_scoped_release>())
    // .def_readonly("v_data", &MyClass::v_data, byref)
    // .def_readonly("v_gamma", &MyClass::v_gamma, byref)
    // ;

    py::class_<Nexullance_IT_interface>(m, "Nexullance_IT_interface")
   .def(py::init<int, Eigen::MatrixX2i, Eigen::MatrixXf, bool>())
   .def("run", &Nexullance_IT_interface::run)
   .def("get_max_link_load", &Nexullance_IT_interface::get_max_link_load)
   .def("get_routing_table", &Nexullance_IT_interface::get_routing_table)
   .def("get_num_attempts_step_2", &Nexullance_IT_interface::get_num_attempts_step_2)
   .def("set_parameters", &Nexullance_IT_interface::set_parameters)
    ;
}
