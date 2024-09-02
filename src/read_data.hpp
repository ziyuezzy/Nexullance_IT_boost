#ifndef READ_DATA_HPP
#define READ_DATA_HPP

#include "definitions.hpp"
#include <Eigen/Dense>

Graph read_graphml(std::string graphmlFile, bool print);
Graph read_graph_from_arcs(int V, Eigen::MatrixX2i arcs, bool print);

void read_matrix(std::string filename, bool print, float** result_matrix, size_t dim);

#endif