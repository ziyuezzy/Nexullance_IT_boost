#ifndef LIB_NEXULLANCE_IT_HPP
#define LIB_NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include "read_data.hpp"
#include <Eigen/Dense>


// TODO: a path can be simlified as "Vertex*"

std::tuple<double, float> run_Nexullance_IT_with_paths(std::string input_graph_path, std::string input_matrix_path, int num_step_1);

std::tuple<double, float> run_Nexullance_IT(int V, Eigen::MatrixX2i arcs, Eigen::MatrixXf input_matrix, int num_step_1);

#endif