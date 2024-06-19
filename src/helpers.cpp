#include "definitions.hpp"
#include "Nexullance_IT.hpp"
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <queue>
#include <iostream>
#include <chrono>

void compute_all_shortest_paths_all_s_d(const Graph &G, std::vector<std::vector<Vertex>>** all_s_d_paths,
                                        const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap){

        // iterate over all sources in the boost graph G
        for(size_t src = 0; src < num_vertices(G); src++){
            // dijkstra_shortest_paths_no_color_map(G, src, weight_map(weightmap));
            compute_all_shortest_paths_single_source(G, src,  all_s_d_paths[src], weightmap);
        }

}

void compute_all_shortest_paths_single_source(const Graph &G, Vertex s, std::vector<std::vector<Vertex>>* all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap) {
    
    float dist[num_vertices(G)];
    dijkstra_shortest_paths_no_color_map(G, s, distance_map(boost::make_iterator_property_map(dist, get(boost::vertex_index, G))).weight_map(weightmap));
    
    // Initialize paths from the source vertex
    all_paths[s].push_back({s});
    
    std::queue<Vertex> q;
    q.push(s);
    
    while (!q.empty()) {
        Vertex u = q.front();
        q.pop();

        OutEdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(u, G); ei != ei_end; ++ei) {
            Edge e = *ei;
            Vertex v = boost::target(e, G);

            if (dist[u] + weightmap[e] == dist[v]) {
                if (all_paths[v].empty()) {
                    q.push(v);
                }
                for (const auto &path : all_paths[u]) {
                    std::vector<Vertex> new_path = path;
                    new_path.push_back(v);
                    all_paths[v].push_back(new_path);
                }
            }
        }
    }
}


void compute_all_shortest_paths_single_s_d(const Graph &G, Vertex s, Vertex d, std::vector<std::vector<Vertex>> &result_all_paths,
                                const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap) {
    std::vector<std::vector<Vertex>> all_paths[num_vertices(G)];
    
    assert(s != d);
    // make sure s and d are valid vertices in the graph
    assert(s < num_vertices(G) && d < num_vertices(G));

    float dist[num_vertices(G)];
    dijkstra_shortest_paths_no_color_map(G, s, distance_map(boost::make_iterator_property_map(dist, get(boost::vertex_index, G))).weight_map(weightmap));
    
    // Initialize paths from the source vertex
    all_paths[s].push_back({s});
    std::queue<Vertex> q;
    q.push(s);
    while (!q.empty()) {
        Vertex u = q.front();
        q.pop();

        OutEdgeIterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(u, G); ei != ei_end; ++ei) {
            Edge e = *ei;
            Vertex v = boost::target(e, G);

            if (dist[u] + weightmap[e] == dist[v]) {
                if (all_paths[v].empty()) {
                    q.push(v);
                }
                for (const auto &path : all_paths[u]) {
                    std::vector<Vertex> new_path = path;
                    new_path.push_back(v);
                    all_paths[v].push_back(new_path);
                }
            }
        }
        if (u == d){
            // properly copy the data from all_paths[u] to the result_all_paths
            result_all_paths = all_paths[u]; // This should be a deep copy
            return;
        }
    }
    // if we reach here, then there is no path from s to d
    // raise an error
    std::cout << "No path found from " << s << " to " << d << std::endl;
    assert(false);
}

std::tuple<double, float> run_Nexullance_IT(std::string input_graph_path, std::string input_matrix_path, int num_step_1){
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


    NexullanceIT nexu_it = NexullanceIT(G, const_cast<const float**>(matrix), Cap_link, debug);
    
    auto start = std::chrono::high_resolution_clock::now();
    nexu_it.optimize(num_step_1, 1.0, 1.0, 6, 0.1, 7.0, num_routers);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // Deallocate the memory for matrix
    for (int i = 0; i < num_routers; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;

    return std::make_tuple(elapsed.count(), nexu_it.result_max_loads_step_2.back());
}
