#ifndef MD_NEXULLANCE_IT_HPP
#define MD_NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>

using path_id = size_t;

class MD_Nexullance_IT{
    // _nx_graph: nxGraph, _M_R: np.ndarray, _Cap_remote: float, _verbose:bool=False
    public:
        MD_Nexullance_IT(Graph& _input_graph, const std::vector<float**> _M_Rs, const std::vector<float> _M_R_weights, const float _Cap_remote, const bool _verbose=false);
        ~MD_Nexullance_IT();

        void get_finial_max_load();
        void step_1(float _alpha, float _beta);
        bool step_2(float _alpha, float _beta, float step, float threshold=0.001, int min_attempts=50, int max_attempts=100000);
        void optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, float alpha_step_2, float beta_step_2, 
                        int method_2_min_attempts, int method_2_threshold, int method_2_max_attempts);

        std::list<float> result_max_loads_step_1;
        std::list<float> result_max_loads_step_2;
        size_t num_attempts_step_2 = 0;

        result_routing_table get_routing_table();
        float get_average_path_length();
        inline std::vector<float> get_max_load_vec(){
            return max_load_vec;
        }
    private:
        Graph G;
        size_t num_edges;
        size_t num_vertices;
        size_t M; // number of demand matrices as input
        const std::vector<float**> M_Rs;
        const std::vector<float> M_R_weights;
        std::vector<float> max_load_vec; // max load under each demand matrix
        const float Cap_remote;
        const bool verbose;
        path_id next_path_id;
        std::unordered_map<path_id, std::vector<Vertex>> path_id_to_path; // using vector here, because the shortest-path algorithm return vector<Vertex> as a path
        std::unordered_map<path_id, float>** routing_tables; // a 2D-array of map, first index corresponds to source router id, second index the destination router id.
        std::vector<float**> link_load_vec;
        std::vector<path_id>** link_path_ids;
        float final_max_load;
        
        // boost::unordered_map<Edge, float> link_load; // TODO?: probably not necessary to use a map. So in case of bottleneck of std::map operation, we can consider to optimize this?(2D s-d array with sparse values)
        // boost::unordered_map<Edge, std::vector<path_id>> link_path_ids; // using vector here, because we may need to dynamically remove or append elements.   

        std::vector<std::vector<Vertex>>** all_paths_all_s_d;
        property_map< Graph, edge_weight_t >::type weightmap;


};
#endif // GRAPH_DEFINITIONS_HPP
