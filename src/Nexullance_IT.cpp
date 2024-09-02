#include "Nexullance_IT.hpp"
#include <boost/graph/adjacency_list.hpp>
// #include <boost/property_map/property_map.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>

NexullanceIT::NexullanceIT(Graph& _input_graph, const float** _M_R, 
    const float _Cap_remote, const bool _verbose): G(_input_graph), 
    M_R(_M_R), Cap_remote(_Cap_remote), verbose(_verbose), final_max_load(NULL) {

    num_edges = boost::num_edges(G);
    num_vertices = boost::num_vertices(G);
    next_path_id = 0;
    routing_tables = new std::unordered_map<path_id, float>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        routing_tables[i] = new std::unordered_map<path_id, float>[num_vertices];
    }

    link_load = new float*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        link_load[i] = new float[num_vertices];
    }
    // assign 0.0 to all link_load initially
    for (int i = 0; i < num_vertices; i++) { // order of loops??
        for (int j = 0; j < num_vertices; j++) {
            link_load[i][j] = 0.0;
        }
    }

    link_path_ids = new std::vector<path_id>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        link_path_ids[i] = new std::vector<path_id>[num_vertices];
    }


    // initialize data structures
    all_paths_all_s_d = new std::vector<std::vector<Vertex>>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        all_paths_all_s_d[i] = new std::vector<std::vector<Vertex>>[num_vertices];
    }
    weightmap = get(edge_weight, G);

    // EdgeIterator ei, ei_end;
    // for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++) {
    //     link_load[*ei] = 0;
    //     link_path_ids[*ei] = {};
    // }
}

NexullanceIT::~NexullanceIT() {
    for (int i = 0; i < num_vertices; i++) {
        delete[] routing_tables[i];
        delete[] link_load[i];
        delete[] link_path_ids[i];
        delete[] all_paths_all_s_d[i];
        // delete[] M_R[i];
    }
    // delete[] M_R;
    delete[] routing_tables;
    delete[] link_load;
    delete[] link_path_ids;
    delete[] all_paths_all_s_d;
}

void NexullanceIT::step_1(float _alpha, float _beta) {
    // first calculate the paths

    
    compute_all_shortest_paths_all_s_d(G, all_paths_all_s_d, weightmap);

    if(verbose)
        std::cout<<"step 1: computed all shortest paths all s d"<<std::endl;

    // then set all loads to 0.0, empty path ids, empty path_id_to_path
    path_id_to_path.clear();
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++) {
        link_load[(*ei).m_source][(*ei).m_target] = 0;
        link_path_ids[(*ei).m_source][(*ei).m_target].clear();
    }
    next_path_id = 0;

    if(verbose)
        std::cout<<"step 1: cleared link load and path ids"<<std::endl;
    
    // clear and update routing table, update link load, update path ids
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            routing_tables[s][d].clear();
            if(s==d){
                continue;
            }
            std::vector<std::vector<Vertex>> paths = all_paths_all_s_d[s][d];
            if(verbose)
                std::cout<<"step 1: cleared routing table for s = "<<s<<" d = "<<d<<std::endl;

            float ECMP_weight = 1.0/paths.size();
            for (std::vector<Vertex> path: paths){
                path_id path_id = next_path_id++;
                path_id_to_path[path_id] = path;
                routing_tables[s][d][path_id] = ECMP_weight;
                for (int l = 0; l < path.size() - 1; l++) {
                    Vertex u = path[l];
                    Vertex v = path[l+1];
                    Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                    link_path_ids[e.m_source][e.m_target].push_back(path_id);
                    link_load[e.m_source][e.m_target] += ECMP_weight*M_R[s][d]/Cap_remote; // TODO: to further optimize, divide Cap_remote on the M_R at the beginning?
                    boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta)); // TODO: alternatively, update the weights outside this loop, may lead to speedup
                }
            }
        }
    }

    if(verbose)
        std::cout<<"step 1: updated routing table, update link load, path ids, and weights"<<std::endl;


    // find the max value of link_load (float** , of size num_vertices*num_vertices)
    float max_load = 0.0;
    for (int i = 0; i < num_vertices; i++) { // order of loops??
        for (int j = 0; j < num_vertices; j++) {
            if (link_load[i][j] > max_load) {
                max_load = link_load[i][j];
            }
        }
    }
    if(verbose)
        std::cout<<"step 1: found max link load"<<std::endl;

    // float max_load = 0;
    // for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
    //     if (link_load[*ei] > max_load) {
    //         max_load = link_load[*ei];
    //     }
    // }
    result_max_loads_step_1.push_back(max_load);
    if(verbose)
        std::cout<<"result from step1: " << result_max_loads_step_1.back() << std::endl;
    return;
}

bool NexullanceIT::step_2(float _alpha, float _beta, float step, float threshold, int min_attempts, int max_attempts) {
    
    property_map< Graph, edge_weight_t >::type weightmap = get(edge_weight, G);
    auto rng = std::default_random_engine {};
    int attempts = 0;
    std::list<float> max_loads;
    float max_load = 0;


    while (attempts < max_attempts) {

        float max_load = link_load[0][0];
        for (int i = 0; i < num_vertices; i++) { // order of loops??
            for (int j = 0; j < num_vertices; j++) {
                if (link_load[i][j] > max_load) {
                    max_load = link_load[i][j];
                }
            }
        }

        // // first find the max value of link_load
        // auto maxElement = std::max_element(link_load.begin(), link_load.end(),
        //     [](const std::pair<Edge, float>& a, const std::pair<Edge, float>& b) {
        //         return a.second < b.second;
        //     });
        // max_load = maxElement->second;
        max_loads.push_back(max_load);

        if((attempts > min_attempts) && ((std::accumulate(std::prev(max_loads.end(), min_attempts), max_loads.end(), 0.0f) / (float)min_attempts -max_load)<threshold)){
            if (verbose){
                std::cout<<"step 2: low progress, terminating for step = "<< step <<std::endl;
                std::cout<<"step 2: found max link load" << max_load <<std::endl;
            }
            result_max_loads_step_2.push_back(max_load);
            return true;
        }
        bool success_attempt = false;

    for (int i = 0; i < num_vertices; i++) { // order of loops??
        for (int j = 0; j < num_vertices; j++) {
        // for (auto iter = link_load.begin(); iter != link_load.end(); ++iter) {
            // main function body start from here:
            // if (iter->second < max_load){
            //     continue;
            // }
            
            // Edge e = iter->first;

            if (link_load[i][j] < max_load){
                continue;
            }
            std::vector<path_id> path_ids = link_path_ids[i][j];

            std::multimap<float, path_id, std::greater<float>> sorted_path_ids;
            for (auto old_path_id: path_ids) {
                std::vector<Vertex> old_path = path_id_to_path[old_path_id];
                int src = old_path.front();
                int dst = old_path.back();
                float contribution = routing_tables[src][dst][old_path_id]*M_R[src][dst]/Cap_remote;
                sorted_path_ids.insert(std::make_pair(contribution, old_path_id));
            }

            for (auto item: sorted_path_ids) {
                path_id old_path_id = item.second;
                std::vector<Vertex> old_path = path_id_to_path[old_path_id];
                Vertex src = old_path.front();
                Vertex dst = old_path.back();


                if(verbose){
                    std::cout<<"step 2: starting with old path: " ;
                    for(auto v: old_path)
                        std::cout<<v<<" ";
                    std::cout<<std::endl;
                }
                        
                std::vector<std::vector<Vertex>> all_paths;
                compute_all_shortest_paths_single_s_d(G, src, dst, all_paths, weightmap);

                if(verbose)
                    std::cout<<"step 2: found " << all_paths.size() << " new paths for s = "<<src<<" d = "<<dst<<std::endl;

                std::shuffle(std::begin(all_paths), std::end(all_paths), rng);

                for(std::vector<Vertex> new_path: all_paths){
                    if(new_path != old_path){

                        // // check whether or not there is a max loaded link in the new path
                        // float new_path_max_load = 0.0;
                        // for (int l = 0; l < new_path.size() - 1; l++) {
                        //     if (link_load[new_path[l]][new_path[l+1]] > new_path_max_load) {
                        //         new_path_max_load = link_load[new_path[l]][new_path[l+1]];
                        //     }
                        // }
                        // if(new_path_max_load == max_load){
                        //     continue;
                        // }

                        success_attempt = true;
                        attempts++;

                        if(verbose){
                            std::cout<<"step 2: starting with new path: " ;
                            for(auto v: new_path)
                                std::cout<<v<<" ";
                            std::cout<<std::endl;
                        }
                        
                        //update the paths    
                        float delta_weigth = NULL;
                        std::unordered_map<path_id, float>& current_routing_table = routing_tables[src][dst];

                        auto iter = current_routing_table.find(old_path_id);
                        assert(iter != current_routing_table.end());
                        float old_path_weight = iter->second;

                        bool new_path_found = false;
                        path_id new_path_id = NULL;
                        for (auto it = current_routing_table.begin(); it != current_routing_table.end(); ++it) {
                            if (path_id_to_path[it->first] == new_path) {
                                new_path_found = true;
                                new_path_id = it->first;
                                float prev_path_weight = it->second;
                                delta_weigth = std::min(step, std::min(old_path_weight, 1 - it->second)); 
                                // delta_weigth = std::min(std::min(step, std::min(old_path_weight, 1 - it->second)), Cap_remote*(max_load-new_path_max_load)/M_R[src][dst]); 
                                current_routing_table.erase(it);
                                current_routing_table[new_path_id] = prev_path_weight+delta_weigth;
                                // it->second += delta_weigth;
                                break;
                            }
                        }
                        if(!new_path_found){
                            new_path_id = next_path_id++;
                            path_id_to_path[new_path_id] = new_path;
                            delta_weigth = std::min(step, old_path_weight);
                            // delta_weigth = std::min(std::min(step, old_path_weight), Cap_remote*(max_load-new_path_max_load)/M_R[src][dst]);
                            current_routing_table[new_path_id] = delta_weigth;
                        }

                        // iterate over the new path and update the link load and path ids
                        for (int l = 0; l < new_path.size() - 1; l++) {
                            Vertex u = new_path[l];
                            Vertex v = new_path[l+1];
                            Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                            link_load[e.m_source][e.m_target] += delta_weigth*M_R[src][dst]/Cap_remote; // TODO: to further optimize, divide Cap_remote on the M_R at the beginning?
                            // // compare with max_load, update if this is higher
                            // if (link_load[e.m_source][e.m_target] > (max_load + 0.001)) {
                            //     assert(false &&"should not happen");
                            //     max_load = link_load[e.m_source][e.m_target];
                            // }
                            
                            boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target] ,_beta));
                            if(!new_path_found)
                                link_path_ids[e.m_source][e.m_target].push_back(new_path_id);
                        }

                        current_routing_table[old_path_id] -= delta_weigth;
                        // iterate over the old path and update the link load and path ids
                        for (int l = 0; l < old_path.size() - 1; l++) {
                            Vertex u = old_path[l];
                            Vertex v = old_path[l+1];
                            Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                            link_load[e.m_source][e.m_target] -= delta_weigth*M_R[src][dst]/Cap_remote; // TODO: to further optimize, divide Cap_remote on the M_R at the beginning?
                            boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta));
                        
                            if (current_routing_table[old_path_id] < 0.00001) {
                                auto iter = std::find(link_path_ids[e.m_source][e.m_target].begin(), link_path_ids[e.m_source][e.m_target].end(), old_path_id);
                                if(iter != link_path_ids[e.m_source][e.m_target].end()){
                                    link_path_ids[e.m_source][e.m_target].erase(iter);
                                }
                            }
                        }
                        if (current_routing_table[old_path_id] < 0.00001) {
                            current_routing_table.erase(old_path_id);
                            path_id_to_path.erase(old_path_id);
                        }
                        
    
                        break;
                    }
                }

                if(success_attempt ){
                    break;
                }else{
                    continue;
                }

            }

            if(success_attempt){
                break;
            }else{
                continue;
            }
        }
    }

        if(!success_attempt){
            if(verbose)
                std::cout<<"step 2: found max link load" << max_load <<std::endl;
            std::cout<<"step 2: no progress, terminating after" << attempts << " attempts"<<std::endl;
            result_max_loads_step_2.push_back(max_load);
            return false;
        }
    }

    if(verbose)
        std::cout<<"step 2: found max link load" << max_load <<std::endl;
    std::cout<<"step 2: max number of attemtps reached with step = "<<step<<", threshold = "<<threshold<<", min_attempts = "<<min_attempts<<", max_attempts = "<<max_attempts<<std::endl;
    result_max_loads_step_2.push_back(max_load);
    return true;
}

void NexullanceIT::optimize(int num_step_1, float alpha_step_1, float beta_step_1, int max_num_step_2, float alpha_step_2, float beta_step_2, int method_2_min_attempts){
    assert(num_step_1 >= 1 and "num_step_1 should be a positive integer greater than 1.");
    for (int i = 0; i < num_step_1; i++) {
        step_1(alpha_step_1, beta_step_1);
    }
    float step = 0.5;
    for (int i = 0; i < max_num_step_2; i++) {
        if(step_2(alpha_step_2, beta_step_2, step, 0.001, method_2_min_attempts)){
            step *= 0.5;
            continue;
        }else{
            break;
        }
    }
}
