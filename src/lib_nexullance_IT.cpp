#include "lib_nexullance_IT.hpp"

#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// ===========class Nexullance_IT_interface==========
// The interface object corresponds to a graph, so it can run (or profile) IT or MD_IT for different input demand matrices
Nexullance_IT_interface::Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, 
    const float Cap_core,const float Cap_access, bool debug): 
    _V(V), _Cap_core(Cap_core), _Cap_access(Cap_access), _debug(debug) {

    G = read_graph_from_arcs(V, arcs, false);
    int num_routers=boost::num_vertices(G);
    assert(_V == num_routers);
}

Nexullance_IT_interface::~Nexullance_IT_interface(){
}

IT_outputs Nexullance_IT_interface::run_IT(Eigen::MatrixXf M_EP, const int EPR){
    // converting the Eigen matrix to a float** matrix (for performance)
    int num_EPs = _V*EPR;
    float** matrix = new float*[num_EPs];
    for (int i = 0; i < num_EPs; i++) {
        matrix[i] = new float[num_EPs];
    }
    assert(M_EP.rows() == num_EPs && M_EP.cols() == num_EPs);
    for (int i = 0; i < num_EPs; ++i) {
        for (int j = 0; j < num_EPs; ++j) {
            matrix[i][j] = M_EP(i, j);
        }
    }
    //===========
    Nexullance_IT nexu_it = Nexullance_IT(G, const_cast<const float**>(matrix), EPR, _Cap_core, _Cap_access, _debug);
    
    // delete the float** matrix
    for (int i = 0; i < num_EPs; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;

    auto start = std::chrono::high_resolution_clock::now();
    nexu_it.optimize(1, 1.0, 1.0, _num_steppings+1, _alpha, _beta, _min_attempts, _stepping_threshold, _max_attempts); // TODO: pass all params
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    IT_outputs result = IT_outputs(elapsed.count(), nexu_it.get_max_core_load(), nexu_it.get_phi(), nexu_it.get_routing_table(), nexu_it.num_attempts_step_2); //TODO: implemnet

    return result;
}

MD_IT_outputs Nexullance_IT_interface::run_MD_IT(std::vector<Eigen::MatrixXf> M_EPs_s, std::vector<float> weights, const int EPR){
    int _M = M_EPs_s.size();
    assert(_M == weights.size());
    int num_EPs = _V*EPR;

    // converting the Eigen matrix to a float** matrix
    std::vector<float**> M_EP_matrices;
    for (int m = 0; m < _M; m++) {
        M_EP_matrices.push_back(new float*[num_EPs]);
        for (int j = 0; j < num_EPs; j++) {
            M_EP_matrices[m][j] = new float[num_EPs];
        }
    }

    for (int m = 0; m < _M; m++) {
        assert(M_EPs_s[m].rows() == num_EPs && M_EPs_s[m].cols() == num_EPs);
        for (int i = 0; i < num_EPs; ++i) {
            for (int j = 0; j < num_EPs; ++j) {
                M_EP_matrices[m][i][j] = M_EPs_s[m](i, j);
            }
        }
    }
    //===========
    MD_Nexullance_IT md_nexu_it = MD_Nexullance_IT(G, M_EP_matrices, weights, EPR, 
                                    _Cap_core, _Cap_access, _debug);

    auto start = std::chrono::high_resolution_clock::now();
    md_nexu_it.optimize(1, 1.0, 1.0, _num_steppings+1, _alpha, _beta, _min_attempts, 
                        _stepping_threshold, _max_attempts, _cal_least_margins);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    MD_IT_outputs result = MD_IT_outputs(elapsed.count(), md_nexu_it.get_max_core_load_vec(), md_nexu_it.get_phis(), md_nexu_it.get_routing_table(), md_nexu_it.get_obj());

    return result;
}
// ==================================================
// ============diff_Nexullance_IT_interface==========

diff_Nexullance_IT_interface::diff_Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, 
    const float Cap_core,const float Cap_access, bool online_mode, bool debug): 
    _V(V), _Cap_core(Cap_core), _Cap_access(Cap_access), _debug(debug), _online_mode(online_mode) {

    G = read_graph_from_arcs(V, arcs, false);
    int num_routers=boost::num_vertices(G);
    assert(_V == num_routers);
}

diff_Nexullance_IT_interface::~diff_Nexullance_IT_interface(){
    if (diff_nexu_it)
        delete diff_nexu_it;
}

std::vector<IT_outputs> diff_Nexullance_IT_interface::run_for_batch_matrices(std::vector<Eigen::MatrixXf> M_EPs_s, const int EPR){
    assert(_online_mode == false);
    int _M = M_EPs_s.size();
    this->_initialize(EPR); // initialize the diff_Nexullance_IT object
    std::vector<IT_outputs> results;
    results.clear();
    for (int m = 0; m < _M; m++)
        results.push_back( this->add_next_matrix(M_EPs_s[m]) );
    assert(results.size() == _M && "results size should be equal to the number of matrices");
    return results;

    // int _M = M_EPs_s.size();
    // int num_EPs = _V*EPR;
    // // converting the Eigen matrix to a float** matrix
    // std::vector<float**> M_EP_matrices;
    // for (int m = 0; m < _M; m++) {
    //     M_EP_matrices.push_back(new float*[num_EPs]); // TODO: delete the new float**
    //     for (int j = 0; j < num_EPs; j++) {
    //         M_EP_matrices[m][j] = new float[num_EPs]; // TODO: delete the new float**
    //     }
    // }
    // for (int m = 0; m < _M; m++) {
    //     assert(M_EPs_s[m].rows() == num_EPs && M_EPs_s[m].cols() == num_EPs);
    //     for (int i = 0; i < num_EPs; ++i)
    //         for (int j = 0; j < num_EPs; ++j)
    //             M_EP_matrices[m][i][j] = M_EPs_s[m](i, j);
    // }
    // //===========
    // diff_nexu_it = new diff_Nexullance_IT(G, EPR, _Cap_core, _Cap_access, _debug, _online_mode);
    // diff_nexu_it->run_for_M_EPs_s(M_EP_matrices, _alpha, _beta, _stepping_threshold, _min_step, _min_attempts, _max_attempts);

}

void diff_Nexullance_IT_interface::_initialize(const int EPR){
    num_EPs = _V*EPR;
    diff_nexu_it = new diff_Nexullance_IT(G, EPR, _Cap_core, _Cap_access, _debug);
    // TODO: set params?
}

IT_outputs diff_Nexullance_IT_interface::add_next_matrix(Eigen::MatrixXf M_EPs){
    // converting the Eigen matrix to a float** matrix, assert that the size is correct
    float** matrix = new float*[num_EPs];
    for (int i = 0; i < num_EPs; i++)
        matrix[i] = new float[num_EPs];

    assert(M_EPs.rows() == num_EPs && M_EPs.cols() == num_EPs);
    for (int i = 0; i < num_EPs; ++i)
        for (int j = 0; j < num_EPs; ++j)
            matrix[i][j] = M_EPs(i, j);
    //===========
    IT_outputs result = diff_nexu_it->optimize_for_M_EPs(matrix, _alpha, _beta, _stepping_threshold, _min_step, _min_attempts, _max_attempts);
    delete[] matrix;
    return result;
}

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(Nexullance_IT_cpp, m) {
    m.doc() = "calling Nexullance IT";

    // m.def("run_Nexullance_IT", &run_Nexullance_IT, "test_func_directly_calling_nexullance_IT"); // example of defining a function for python

    py::class_<Nexullance_IT_interface>(m, "Nexullance_IT_interface")
   .def(py::init<int, Eigen::MatrixX2i, const float, const float, bool>(), py::arg("V"), py::arg("arcs"), py::arg("Cap_core") = 10, py::arg("Cap_access") = 10, py::arg("debug") = false)
   .def("set_parameters", &Nexullance_IT_interface::set_parameters)
   .def("run_IT", &Nexullance_IT_interface::run_IT)
   .def("run_MD_IT", &Nexullance_IT_interface::run_MD_IT)
   ;

    py::class_<diff_Nexullance_IT_interface>(m, "diff_Nexullance_IT_interface")
    .def(py::init<int, Eigen::MatrixX2i, const float, const float, bool, bool>(), py::arg("V"), py::arg("arcs"), py::arg("Cap_core") = 10, py::arg("Cap_access") = 10, py::arg("online_mode") = false, py::arg("debug") = false)
    .def("set_parameters", &diff_Nexullance_IT_interface::set_parameters)
    .def("run_for_batch_matrices", &diff_Nexullance_IT_interface::run_for_batch_matrices)
    .def("_initialize", &diff_Nexullance_IT_interface::_initialize)
    .def("add_next_matrix", &diff_Nexullance_IT_interface::add_next_matrix)
    ;

    py::class_<IT_outputs>(m, "IT_outputs")
        .def("get_elapsed_time", &IT_outputs::get_elapsed_time)
        .def("get_max_core_link_load", &IT_outputs::get_max_core_link_load)
        .def("get_phi", &IT_outputs::get_phi)
        .def("get_routing_table", &IT_outputs::get_routing_table)
        .def("get_num_attempts", &IT_outputs::get_num_attempts);

    py::class_<MD_IT_outputs>(m, "MD_IT_outputs")
        .def("get_elapsed_time", &MD_IT_outputs::get_elapsed_time)
        .def("get_max_core_link_load", &MD_IT_outputs::get_max_core_link_loads)
        .def("get_phis", &MD_IT_outputs::get_phis)
        .def("get_routing_table", &MD_IT_outputs::get_routing_table)
        .def("get_obj", &MD_IT_outputs::get_obj);
        // .def("get_weighted_sum_phi", &MD_IT_outputs::get_weighted_sum_phi)

}
