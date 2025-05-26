#include "MD_Nexullance_IT.hpp"
#include "graph_utility.hpp"
#include <boost/graph/adjacency_list.hpp>
// #include <boost/property_map/property_map.hpp>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <random>

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

MD_Nexullance_IT::MD_Nexullance_IT(Graph& _input_graph, std::vector<float**> _M_EPs_s, const std::vector<float> _M_weights, 
                        const int _EPR, const float _Cap_core, const float _Cap_access, const bool _verbose): 
                        G(_input_graph), Cap_core(_Cap_core), Cap_access(_Cap_access), verbose(_verbose), EPR(_EPR), M_weights(_M_weights){

    num_edges = boost::num_edges(G);
    num_vertices = boost::num_vertices(G);

    M = _M_EPs_s.size();
    assert(M == M_weights.size());

    // Ensure the weights sum to 1
    float sum = std::accumulate(M_weights.begin(), M_weights.end(), 0.0f);
    assert(std::abs(sum - 1.0f) < 1e-6 && "Error: Weights do not sum to 1!");
    // Create a vector of cumulative sums
    cumulative_weights.resize(M_weights.size());
    std::partial_sum(M_weights.begin(), M_weights.end(), cumulative_weights.begin());


    M_Rs.resize(M);
    total_flows_per_EP.resize(M);
    max_access_loads.resize(M);

    for (size_t i = 0; i < M; i++)
    {
        std::pair<float, float> temp_result=procress_M_EPs(const_cast<const float**>(_M_EPs_s[i]), num_vertices, EPR, &(M_Rs[i]));
        // Calculate the max_access_link_load and total flow for later usage
        total_flows_per_EP[i] = temp_result.second/(num_vertices*EPR);
        max_access_loads[i] = temp_result.first/Cap_access;
    }
    
    max_core_loads.resize(M, 0.0f);
    phis.resize(M, 0.0f);
    // Obj_map_S.clear();
    // weighted_sum_phi = 0.0f;
    Objective_func_inverse = 0.0f;

    max_core_load_links.resize(M);

    next_path_id = 0;
    routing_tables = new std::unordered_map<path_id, float>*[num_vertices];
    for (int i = 0; i < num_vertices; i++) {
        routing_tables[i] = new std::unordered_map<path_id, float>[num_vertices];
    }

    link_load_vec.clear();
    link_load_vec.reserve(M);
    for (int m = 0; m < M; m++) {
        link_load_vec.push_back(new float*[num_vertices]);
        for (int j = 0; j < num_vertices; j++) {
            link_load_vec[m][j] = new float[num_vertices]{0.0f};
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
}

MD_Nexullance_IT::~MD_Nexullance_IT() {
    for (int i = 0; i < num_vertices; i++) {
        delete[] routing_tables[i];
        delete[] link_path_ids[i];
        delete[] all_paths_all_s_d[i];
    }
    for (int m = 0; m < M; m++) {   
        for (int i = 0; i < num_vertices; i++) {
            delete[] link_load_vec[m][i];
            delete[] M_Rs[m][i];
        }
        delete[] link_load_vec[m];
        delete[] M_Rs[m];
    }
    
    delete[] routing_tables;
    delete[] link_path_ids;
    delete[] all_paths_all_s_d;
}

void MD_Nexullance_IT::step_1(float _alpha, float _beta) {
    // first calculate the paths

    compute_all_shortest_paths_all_s_d(G, all_paths_all_s_d, weightmap);
    // if(verbose)
    //     std::cout<<"step 1: computed all shortest paths all s d"<<std::endl;

    // then set all loads to 0.0, empty path ids, empty path_id_to_path
    path_id_to_path.clear();
    EdgeIterator ei, ei_end;
    for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++) {
        for (int m = 0; m < M; m++) {
            link_load_vec[m][(*ei).m_source][(*ei).m_target] = 0.0;
        }
        link_path_ids[(*ei).m_source][(*ei).m_target].clear();
    }
    next_path_id = 0;

    // if(verbose)
    //     std::cout<<"step 1: cleared link load and path ids"<<std::endl;
    
    // clear and update routing table, update link load, update path ids
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            routing_tables[s][d].clear();
            if(s==d){
                continue;
            }
            std::vector<std::vector<Vertex>> paths = all_paths_all_s_d[s][d];
            // if(verbose)
            //     std::cout<<"step 1: cleared routing table for s = "<<s<<" d = "<<d<<std::endl;

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
                    for (int m = 0; m < M; m++) {
                        link_load_vec[m][e.m_source][e.m_target] += ECMP_weight*M_Rs[m][s][d]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                    }
                    // boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta)); // TODO: alternatively, update the weights outside this loop, may lead to speedup
                }
            }
        }
    }

    // if(verbose)
    //     std::cout<<"step 1: updated routing table, update link load, path ids, and weights"<<std::endl;

    // find the max link load of each m, and then calculate the weighted average
    for (int m = 0; m < M; m++) {
        max_core_loads[m] = 0.0;

        EdgeIterator ei, ei_end;
        for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
            size_t i = (*ei).m_source;
            size_t j = (*ei).m_target;
            if (link_load_vec[m][i][j] > max_core_loads[m]) {
                max_core_loads[m] = link_load_vec[m][i][j];
            }
        }
        phis[m] = total_flows_per_EP[m]/std::max(max_access_loads[m], max_core_loads[m]);
        // weighted_sum_phi += phis[m]*M_weights[m];
        Objective_func_inverse += M_weights[m]/phis[m];
    }

    if(verbose)
        std::cout<<"result objective function from step1: " << 1/Objective_func_inverse << std::endl;
    return;
}

bool MD_Nexullance_IT::step_2(float _alpha, float _beta, float step, float threshold, int min_attempts, int max_attempts, bool cal_least_margins) {
    
    if (verbose)
        std::cout<<"step 2: starting for step = "<< step <<std::endl;
        
    auto rng = std::default_random_engine {};
    int attempts = 0;
    std::list<float> Objective_func_inverse_hist;

    // TODO: now it only keep the most recent attempt in memory, maybe also include further attempts in the history?
    long last_decreased_path_id = -1;
    long last_increased_path_id = -1;

    while (attempts < max_attempts) {
        Objective_func_inverse = 0.0;
        // find the max link load of each m, and then calculate the weighted average
        for (int m = 0; m < M; m++) {
            max_core_load_links[m].clear();
            max_core_loads[m] = 0.0;


            EdgeIterator ei, ei_end;
            for (tie(ei, ei_end) = boost::edges(G); ei != ei_end; ei++) {
                size_t i = (*ei).m_source;
                size_t j = (*ei).m_target;
                if (link_load_vec[m][i][j] > max_core_loads[m]) {
                    max_core_loads[m] = link_load_vec[m][i][j];
                    max_core_load_links[m].clear();
                    max_core_load_links[m].push_back(std::make_pair(i,j));
                } else if (link_load_vec[m][i][j] == max_core_loads[m]) {
                    max_core_load_links[m].push_back(std::make_pair(i,j));
                }
            }
            // float Gamma_m = max_core_loads[m]*M_R_weights[m]; // the products of max_core_loads[m] and M_R_weights[m]
            // weighted_sum_phi += Gamma_m;  
            phis[m] = total_flows_per_EP[m]/std::max(max_access_loads[m], max_core_loads[m]);
            float weighted_contribution = M_weights[m]/phis[m];
            // Obj_map_S.insert(std::make_pair(weighted_contribution, m));
            Objective_func_inverse += weighted_contribution;
        }
        Objective_func_inverse_hist.push_back(Objective_func_inverse);
        // assert(n != -1);

        // now determine n (which demand matrix to work on)

        // // =================simple round robin=============
        // int n = attempts%M; // round robin starting point
        // int counter = 0;
        // while (counter < M) // loop through all demand matrices at most once
        // {
        //     if(max_access_loads[n] < max_core_loads[n]){ 
        //         break; // if the core link load is larger, we choose this demand matrix and therefore break the loop
        //     }
        //     // if the access link load becomes larger, skip this demand matrix, try the next one
        //     if(verbose){
        //         std::cout<<"access load becomes bottleneck, skipping demand matrix "<< n << std::endl;
        //     }  
        //     n = (n+1)%M;
        //     counter+=1;
        // }
        // // if all m satisfies "max_access_loads[m] >= max_core_loads[m]", terminate.
        // if (counter == M){
        //     if(verbose){
        //         std::cout<<"access load becomes bottleneck in all demand matrices, terminating with objective function = " << 1/Objective_func_inverse << std::endl;
        //     }  
        //     return false;
        // }
        // // =================================================
        // // ==============iterate through Obj_map_S==========
        // uint n;
        // auto it = Obj_map_S.begin();
        // while (it != Obj_map_S.end()) // loop through all demand matrices at most once
        // {
        //     n = it->second;
        //     if(max_access_loads[n] < max_core_loads[n]){ 
        //         break; // if the core link load is larger, we choose this demand matrix and therefore break the loop
        //     }else{
        //         // if the access link load becomes larger, skip this demand matrix, try the next one
        //         if(verbose){
        //             std::cout<<"access load becomes bottleneck, skipping demand matrix "<< n << std::endl;
        //         }  
        //         n = (n+1)%M;
        //         it++;
        //     }
        // }
        // // if all m satisfies "max_access_loads[m] >= max_core_loads[m]", terminate.
        // if (it == Obj_map_S.end()){
        //     if(verbose){
        //         std::cout<<"access load becomes bottleneck in all demand matrices, terminating with objective function = " << 1/Objective_func_inverse  << std::endl;
        //     }  
        //     return false;
        // }
        // // =================================================
        // =================randomization with weights=============
        // check whether or not all demand matrix satisfy "max_access_loads[n] >= max_core_loads[n]"
        bool should_exit = true;
        for (size_t m = 0; m < M; m++){
            if (max_access_loads[m] < max_core_loads[m]){
                should_exit = false;
            }
        }
        if(should_exit){
            std::cout<<"access load becomes bottleneck in all demand matrices, terminating with objective function = " << 1/Objective_func_inverse  << std::endl;
            return false;
        }
        // throw dice(s) to determine n.
        uint n;
        // Create a distribution that generates numbers between 0 and 1
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        do
        {
            float random_float = dist(rng);
            // Determine which range the random number falls into
            for (size_t i = 0; i < cumulative_weights.size(); ++i) {
                if (random_float < cumulative_weights[i]) {
                    n=i;
                    break;
                }
            }
        } while (max_access_loads[n] >= max_core_loads[n]);
        // =================================================

        if(verbose){
            std::cout<<"step2, it="<< attempts << ", Objective function = " << 1/Objective_func_inverse << std::endl;
            std::cout<<"attempting on n = "<< n << std::endl;
        }        

        if((attempts > min_attempts) && ((std::accumulate(std::prev(Objective_func_inverse_hist.end(), min_attempts/2), Objective_func_inverse_hist.end(), 0.0f)/((float)min_attempts/2) - Objective_func_inverse)<threshold)){
            if (verbose){
                std::cout<<"step 2: low progress, terminating for step = "<< step <<std::endl;
                std::cout<<"average of previous "<< min_attempts/2 <<" steps = "<< 
                    std::accumulate(std::prev(Objective_func_inverse_hist.end(), min_attempts/2), Objective_func_inverse_hist.end(), 0.0f)/((float)min_attempts/2) <<std::endl;
                // std::cout<<"step 2: found max link load" << max_core_loads[n] << " for demand matrix " << n <<std::endl;
            }
            num_attempts_step_2 += attempts;
            return true;
        }

        bool success_attempt = false;

        std::vector<float> new_path_max_loads(M, 0.0);
        for(auto link: max_core_load_links[n]){ // re-route paths regarding the n^th demand matrix
                int i = link.first;
                int j = link.second;

                #ifdef DEBUG
                assert(link_load_vec[n][i][j] == max_core_loads[n]);
                #endif

                std::vector<path_id> path_ids = link_path_ids[i][j];

                std::multimap<float, path_id, std::greater<float>> sorted_path_ids;
                for (auto old_path_id: path_ids) {
                    std::vector<Vertex> old_path = path_id_to_path[old_path_id];
                    int src = old_path.front();
                    int dst = old_path.back();
                    float contribution = routing_tables[src][dst][old_path_id]*M_Rs[n][src][dst]/Cap_core;
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
                    // set edge weights in the graph according to the n^th demand matrix
                    EdgeIterator ei, ei_end;
                    for (tie(ei, ei_end) = boost::edges(G); ei!= ei_end; ei++)
                        boost::put(boost::edge_weight, G, *ei, _alpha + pow(link_load_vec[n][(*ei).m_source][(*ei).m_target],_beta));
                    compute_all_shortest_paths_single_s_d(G, src, dst, all_paths, weightmap);

                    if(verbose)
                        std::cout<<"step 2: found " << all_paths.size() << " new paths for s = "<<src<<" d = "<<dst<<std::endl;

                    std::shuffle(std::begin(all_paths), std::end(all_paths), rng);

                    for(std::vector<Vertex> new_path: all_paths){
                        if(new_path != old_path){

                            // check the max loaded link in the new path for every demand matrix
                            new_path_max_loads.assign(M, 0.0);
                            for (int m = 0; m < M; m++)
                                for (int l = 0; l < new_path.size() - 1; l++)
                                    if (link_load_vec[m][new_path[l]][new_path[l+1]] > new_path_max_loads[m])
                                        new_path_max_loads[m] = link_load_vec[m][new_path[l]][new_path[l+1]];

                            if(new_path_max_loads[n] == max_core_loads[n])  continue;
                            #ifdef DEBUG
                            assert(new_path_max_loads[n] < max_core_loads[n]);
                            #endif
                            std::vector<float> new_path_margins;
                            float least_margin=1.0;
                            if (cal_least_margins){
                                new_path_margins.clear();
                                for (int m = 0; m < M; m++)
                                    if (M_Rs[m][src][dst]>0)
                                        new_path_margins.push_back(Cap_core*(max_core_loads[m]-new_path_max_loads[m])/M_Rs[m][src][dst]);

                                if (new_path_margins.size() > 0)
                                    least_margin = *std::min_element(new_path_margins.begin(), new_path_margins.end());
                                else
                                    least_margin = 1.0;
                            }

                            if(verbose){
                                std::cout<<"step 2: starting with new path: " ;
                                for(auto v: new_path)
                                    std::cout<<v<<" ";
                                std::cout<<std::endl;
                            }
                            if(verbose && cal_least_margins){
                                std::cout<<"least margin = " << least_margin << std::endl;
                            }
                            //update the paths    
                            float delta_weigth = -1.0;
                            std::unordered_map<path_id, float>& current_routing_table = routing_tables[src][dst];

                            auto iter = current_routing_table.find(old_path_id);
                            #ifdef DEBUG
                            assert(iter != current_routing_table.end());
                            #endif
                            float old_path_weight = iter->second;

                            bool new_path_found = false;
                            long new_path_id = -1;
                            float prev_path_weight = -1.0;
                            for (auto it = current_routing_table.begin(); it != current_routing_table.end(); ++it) {
                                if (path_id_to_path[it->first] == new_path) {
                                    new_path_found = true;
                                    new_path_id = it->first;
                                    prev_path_weight = it->second;
                                    break;
                                }
                            }

                            if (new_path_found){
                                #ifdef DEBUG
                                assert(new_path_id != -1);
                                assert(prev_path_weight != -1.0);
                                #endif
                                if ((old_path_id == last_increased_path_id) && (new_path_id == last_decreased_path_id)){
                                    if(verbose){
                                        std::cout<<"step 2: stopped with the new path for avoiding deadlock "<<std::endl;
                                    }
                                    continue;
                                }
                                if (cal_least_margins){
                                    delta_weigth = std::min(std::min(step, std::min(old_path_weight, 1 - prev_path_weight)), least_margin); 
                                }else{
                                    // still calculate the least margin for the n^th demand matrix
                                    delta_weigth = std::min(std::min(step, std::min(old_path_weight, 1 - prev_path_weight)), 
                                                    Cap_core*(max_core_loads[n]-new_path_max_loads[n])/M_Rs[n][src][dst]); 
                                }
                                current_routing_table.erase(new_path_id);
                                current_routing_table.insert(std::make_pair(new_path_id, prev_path_weight+delta_weigth));

                            }else{
                                #ifdef DEBUG
                                assert(prev_path_weight == -1.0);
                                #endif
                            
                                new_path_id = next_path_id++;
                                path_id_to_path[new_path_id] = new_path;
                                if (cal_least_margins){
                                    delta_weigth = std::min(std::min(step, old_path_weight), least_margin); 
                                }else{
                                    delta_weigth = std::min(std::min(step, old_path_weight), 
                                                    Cap_core*(max_core_loads[n]-new_path_max_loads[n])/M_Rs[n][src][dst]); 
                                }
                                current_routing_table.insert(std::make_pair(new_path_id, delta_weigth));                
                                // current_routing_table[new_path_id] = delta_weigth;
                            }

                            success_attempt = true;
                            attempts++;
                            #ifdef DEBUG
                            assert(delta_weigth != -1.0);
                            #endif

                            last_decreased_path_id = old_path_id;
                            last_increased_path_id = new_path_id;

                            // iterate over the new path and update the link load and path ids
                            for (int l = 0; l < new_path.size() - 1; l++) {
                                Vertex u = new_path[l];
                                Vertex v = new_path[l+1];
                                Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                                if(!new_path_found)
                                    link_path_ids[e.m_source][e.m_target].push_back(new_path_id);
                                for (int m = 0; m < M; m++) { // for each demand matrix
                                    link_load_vec[m][e.m_source][e.m_target] += delta_weigth*M_Rs[m][src][dst]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                                }
                            }

                            current_routing_table.erase(old_path_id);
                            current_routing_table.insert(std::make_pair(old_path_id, old_path_weight-delta_weigth));

                            // iterate over the old path and update the link load and path ids
                            for (int l = 0; l < old_path.size() - 1; l++) {
                                Vertex u = old_path[l];
                                Vertex v = old_path[l+1];
                                Edge e = boost::edge(u, v, G).first; // the boost::edge(u, v, G) returns a pair<Edge edge(u,v), bool found>
                                for (int m = 0; m < M; m++) { // for each demand matrix
                                    link_load_vec[m][e.m_source][e.m_target] -= delta_weigth*M_Rs[m][src][dst]/Cap_core; // TODO: to further optimize, divide Cap_core on the M_R at the beginning?
                                    // handle rounding errors:
                                    if (link_load_vec[m][e.m_source][e.m_target] < 0) {
                                        #ifdef DEBUG
                                        assert(link_load_vec[m][e.m_source][e.m_target]>-0.000001);
                                        #endif
                                        // link_load_vec[m][e.m_source][e.m_target] = 0;
                                    }                                    // boost::put(boost::edge_weight, G, e, _alpha + pow(link_load[e.m_source][e.m_target],_beta));
                                }
                                if (current_routing_table[old_path_id] < 0.000001) {
                                    auto iter = std::find(link_path_ids[e.m_source][e.m_target].begin(), link_path_ids[e.m_source][e.m_target].end(), old_path_id);
                                    if(iter != link_path_ids[e.m_source][e.m_target].end()){
                                        link_path_ids[e.m_source][e.m_target].erase(iter);
                                    }
                                }
                            }
                            if (current_routing_table[old_path_id] < 0.000001) {
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
        // }

        if(!success_attempt){
            if(verbose){
                std::cout<<"step 2: reaches Objective function = " << 1/Objective_func_inverse <<std::endl;
                std::cout<<"step 2: no progress, terminating after" << attempts << " attempts"<<std::endl;
                }
            // result_max_loads_step_2.push_back(Objective_func_inverse_hist.back());
            num_attempts_step_2 += attempts;
            return false;
        }
    }

    if(verbose){
        std::cout<<"step 2: reaches Objective function =  " << 1/Objective_func_inverse_hist.back() <<std::endl;
        std::cout<<"step 2: max number of attemtps reached with step = "<<step<<", threshold = "<<threshold<<", min_attempts = "<<min_attempts<<", max_attempts = "<<max_attempts<<std::endl;
        }
    // result_max_loads_step_2.push_back(Objective_func_inverse_hist.back());
    num_attempts_step_2 += attempts;
    return true;
}

void MD_Nexullance_IT::optimize(int num_step_1, float alpha_step_1, float beta_step_1, 
                        int max_num_step_2, float alpha_step_2, float beta_step_2, 
                        int method_2_min_attempts, float method_2_threshold, int method_2_max_attempts, bool cal_least_margins){
    assert(num_step_1 >= 1 and "num_step_1 should be a positive integer greater than 1.");
    for (int i = 0; i < num_step_1; i++) {
        step_1(alpha_step_1, beta_step_1);
    }
    float step = 0.5;
    for (int i = 0; i < max_num_step_2; i++) {
        bool _continues = step_2(alpha_step_2, beta_step_2, step, method_2_threshold, method_2_min_attempts, method_2_max_attempts, cal_least_margins);
        if(_continues){
        // if(step_2(alpha_step_2, beta_step_2, step, 0.001, method_2_min_attempts)){
            step *= 0.5;
            continue;
        }else{
            break;
        }
    }
}

result_routing_table MD_Nexullance_IT::get_routing_table(){
    result_routing_table result = result_routing_table();
    // iterate over "routing_tables" and convert it to "result_routing_table"
    for (int s = 0; s < num_vertices; s++) {
        for (int d = 0; d < num_vertices; d++) {
            if (s==d)
                continue;

            std::vector< std::pair< std::vector<Vertex>,float> > paths = std::vector< std::pair< std::vector<Vertex>,float> >();
            for (auto item: routing_tables[s][d]) {
                std::vector<Vertex> path = path_id_to_path[item.first];
                float weight = item.second;
                paths.push_back(std::make_pair(path, weight));
            }
            result.insert(std::make_pair(std::make_pair(s,d), paths));
        }
    }
    return result;
}

// float MD_Nexullance_IT::get_average_path_length(){
//     float ave = 0;
//     for (int s = 0; s < num_vertices; s++) {
//         for (int d = 0; d < num_vertices; d++) {
//             if (s==d)
//                 continue;
//             for (auto item: routing_tables[s][d]) {
//                 std::vector<Vertex> path = path_id_to_path[item.first];
//                 float weight = item.second;
//                 ave += path.size()*weight;
//             }
//         }
//     }
//     ave = ave/num_vertices/num_vertices;
//     return ave;
// }

