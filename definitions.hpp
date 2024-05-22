#ifndef GRAPH_DEFINITIONS_HPP
#define GRAPH_DEFINITIONS_HPP
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
using namespace boost;

// Define a directed graph with integer vertex index and float edge weights
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
                               boost::property<boost::vertex_index_t, int>,
                               boost::property<boost::edge_weight_t, float>
                              > Graph;

// // Define a directed graph with integer vertex index, float edge weights, and integer edge indices
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
//                                boost::property<boost::vertex_index_t, int>,
//                                boost::property<boost::edge_weight_t, float>,
//                                boost::property<boost::edge_index_t, int>
//                               > Graph;


typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIterator;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIterator;

Graph read_graphml(std::string graphmlFile, bool print);

void read_matrix(std::string filename, bool print, float** result_matrix, size_t dim);

// TODO: a path can be simlified as "Vertex*"
void compute_all_shortest_paths_all_s_d(const Graph &G, std::vector<std::vector<Vertex>>** all_s_d_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void compute_all_shortest_paths_single_source(const Graph &G, Vertex s, std::vector<std::vector<Vertex>>* all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void compute_all_shortest_paths_single_s_d(const Graph &G, Vertex s, Vertex d, std::vector<std::vector<Vertex>> &all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap);

void test_shortest_paths(std::string input_graph_path);

std::tuple<double, float> run_Nexullance_IT(std::string input_graph_path, std::string input_matrix_path, int num_step_1);


#endif // GRAPH_DEFINITIONS_HPP
