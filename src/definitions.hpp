#ifndef GRAPH_DEFINITIONS_HPP
#define GRAPH_DEFINITIONS_HPP
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>

// #include <Eigen/Dense>
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

typedef std::map<std::pair<Vertex,Vertex>, std::vector< std::pair<std::vector<Vertex>,float> > > result_routing_table;

void test_shortest_paths(std::string input_graph_path);


#endif // GRAPH_DEFINITIONS_HPP
