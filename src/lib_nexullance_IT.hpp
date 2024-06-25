#ifndef LIB_NEXULLANCE_IT_HPP
#define LIB_NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include "read_data.hpp"
#include "Nexullance_IT.hpp"
#include "MD_Nexullance_IT.hpp"
#include <Eigen/Dense>


// TODO: a path can be simlified as "Vertex*"

std::tuple<double, float> run_Nexullance_IT_with_paths(std::string input_graph_path, std::string input_matrix_path, int num_step_1);

std::tuple<double, float> run_Nexullance_IT(int V, Eigen::MatrixX2i arcs, Eigen::MatrixXf input_matrix, int num_step_1);

class Nexullance_IT_interface {
    public:
    Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, Eigen::MatrixXf input_matrix, bool debug);
    ~Nexullance_IT_interface();
    double run();
    float get_max_link_load();
    result_routing_table get_routing_table();
    size_t get_num_attempts_step_2();
    void set_parameters(float _alpha, float _beta);

    private:
    Nexullance_IT *nexu_it = nullptr;
    int _V;
    float alpha=0.1;
    float beta=7.0;

};

class MD_Nexullance_IT_interface {
    public:
    MD_Nexullance_IT_interface(int V, Eigen::MatrixX2i arcs, std::vector<Eigen::MatrixXf> MRs,
                               std::vector<float> MR_weights, bool debug);
    ~MD_Nexullance_IT_interface();
    double run();
    float get_weighted_max_link_load();
    std::vector<float> get_max_link_loads();
    result_routing_table get_routing_table();
    size_t get_num_attempts_step_2();
    void set_parameters(float _alpha, float _beta);

    private:
    MD_Nexullance_IT *md_nexu_it = nullptr;
    int _V;
    float alpha=0.1;
    float beta=7.0;
};

#endif