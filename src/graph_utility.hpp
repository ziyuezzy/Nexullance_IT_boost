#ifndef GRAPH_UTILITY_HPP
#define GRAPH_UTILITY_HPP

#include "definitions.hpp"

// TODO: a path can be simlified as "Vertex*"
void compute_all_shortest_paths_all_s_d(const Graph &G, std::vector<std::vector<Vertex>>** all_s_d_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void compute_all_shortest_paths_single_source(const Graph &G, Vertex s, std::vector<std::vector<Vertex>>* all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void compute_all_shortest_paths_single_s_d(const Graph &G, Vertex s, Vertex d, std::vector<std::vector<Vertex>> &all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

#endif