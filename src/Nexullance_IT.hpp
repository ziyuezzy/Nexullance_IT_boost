#ifndef NEXULLANCE_IT_HPP
#define NEXULLANCE_IT_HPP

#include "definitions.hpp"
#include <boost/unordered/unordered_map.hpp>
// #include <list>
#include <unordered_map>

using path_id = size_t;

class Nexullance_IT{
    // _nx_graph: nxGraph, _M_R: np.ndarray, _Cap_remote: float, _verbose:bool=False
    public:
        Nexullance_IT(Graph& _input_graph, const float** _M_R, const float _Cap_remote, const bool _verbose=false);
        ~Nexullance_IT();

        void get_finial_max_load();
        void step_1(float _alpha, float _beta);
        bool step_2(float _alpha, float _beta, float step, float threshold=0.001, int min_attempts=50, int max_attempts=100000);
        void optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, float alpha_step_2, float beta_step_2, int method_2_min_attempts);

        std::list<float> result_max_loads_step_1;
        std::list<float> result_max_loads_step_2;
        float final_max_load;
        size_t num_attempts_step_2 = 0;

        typedef std::map<std::pair<Vertex,Vertex>, std::vector< std::pair<std::vector<Vertex>,float> > > result_routing_table;
        result_routing_table get_routing_table();
    private:
        Graph G;
        size_t num_edges;
        size_t num_vertices;
        const float** M_R;
        const float Cap_remote;
        const bool verbose;
        path_id next_path_id;
        std::unordered_map<path_id, std::vector<Vertex>> path_id_to_path; // using vector here, because the shortest-path algorithm return vector<Vertex> as a path
        std::unordered_map<path_id, float>** routing_tables; // a 2D-array of map, first index corresponds to source router id, second index the destination router id.
        float** link_load;
        std::vector<path_id>** link_path_ids;
        
        // boost::unordered_map<Edge, float> link_load; // TODO?: probably not necessary to use a map. So in case of bottleneck of std::map operation, we can consider to optimize this?(2D s-d array with sparse values)
        // boost::unordered_map<Edge, std::vector<path_id>> link_path_ids; // using vector here, because we may need to dynamically remove or append elements.   

        std::vector<std::vector<Vertex>>** all_paths_all_s_d;
        property_map< Graph, edge_weight_t >::type weightmap;


};
#endif // GRAPH_DEFINITIONS_HPP
