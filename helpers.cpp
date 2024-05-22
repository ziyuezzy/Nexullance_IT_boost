#include "definitions.hpp"
#include <boost/graph/dijkstra_shortest_paths_no_color_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <queue>
#include <iostream>

// void compute_all_predecessors(const Graph &G, Vertex s,
//                               const std::vector<Vertex> &dist,
//                               std::list<Vertex> *all_preds,
//                               const boost::property_map<Graph, boost::edge_weight_t>::type &weightmap) {
//     int num_vertices = boost::num_vertices(G);

//     VertexIterator vi, vi_end;
//     for (boost::tie(vi, vi_end) = boost::vertices(G); vi != vi_end; ++vi) {
//         boost::graph_traits<Graph>::vertex_descriptor v = *vi;
//         if (v == s)
//             continue;
//         InEdgeIterator ei, ei_end;
//         for (boost::tie(ei, ei_end) = boost::in_edges(v, G); ei != ei_end; ++ei) {
//             Edge e = *ei;
//             Vertex u = boost::source(e, G);
//             if (dist[u] + weightmap[e] == dist[v]) {
//                 all_preds[v].push_back(u);
//             }
//         }
//     }
// }


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
