import os
import sys
sys.path.append("/users/ziyzhang/topology-research")
from topologies import HPC_topo
import globals as gl
import numpy as np
import pickle
import csv
from nexullance.ultility import nexullance_exp_container
from nexullance import Nexullance_MP

Cap_core = 10 #GBps
Cap_access = 10 #GBps

csv_cols = ['V', 'D', 'M_name', 'M_weight', 'phi_ECMP_ASP', 'phi_MP_APST4', 'phi_IT', 'phi_MD_IT', 'phi_MD_MP_APST4',
            'phi_MP_APST4_uniform', 'phi_MP_APST4_shift_1', 'phi_MP_APST4_halfshift']


def main():
    config = gl.ddf_configs[0]
    V = config[0]
    D = config[1]
    EPR = (D+1)//2
    topo_name="DDF"
    _network = HPC_topo.HPC_topo.initialize_child_instance(topo_name+"topo", V, D)
    arcs = _network.generate_graph_arcs()

    _network.pre_calculate_ECMP_ASP()
    _network.pre_calculate_APST_n(4)

    ECMP_ASP_phis = []
    Nexullance_MP_APST4_uniform_phis = []
    Nexullance_MP_APST4_shift_1_phis = []
    Nexullance_MP_APST4_halfshift_phis = []

    Nexullance_MP_APST4_phis = []
    Nexullance_IT_phis = []

    # ==============================================================
    M_EPs_names=['uniform', 'shift_1', 'half_shift'] # traffic demand matrices to test on
    M_EPs_s = gl.gen_M_EPs_s(topo_name, V, D, EPR, M_EPs_names, 10.0)
    M_weights = [0.1, 0.6, 0.3]


    # ==============================================================
    Nexullance_MP_APST4_M_s = gl.gen_M_EPs_s(topo_name, V, D, EPR, ['uniform', 'shift_1', 'half_shift'], 10.0)

    # calculate the routing table of Nexullance_MP_APST4 for 'uniform'
    nexu = Nexullance_MP.Nexullance_MP(_network.nx_graph, _network.__getattribute__(f"APST_{4}") ,
                                        Cap_core, Cap_access, V, Nexullance_MP_APST4_M_s[0])
    nexu.init_model()
    _, MP_APST4_RT = nexu.solve()
    for M_EPs in M_EPs_s:
        # apply this routing table on all traffic demand matrices
        core_link_flows, access_link_flows = _network.distribute_M_EPs_on_weighted_paths(MP_APST4_RT, EPR, M_EPs)
        # calculate phi
        Nexullance_MP_APST4_uniform_phis.append(gl.network_total_throughput(M_EPs, max(core_link_flows)/Cap_core, max(access_link_flows)/Cap_access)/(V*EPR))
    # =================================================================

    # calculate the routing table of Nexullance_MP_APST4 for 'shift_1'
    nexu = Nexullance_MP.Nexullance_MP(_network.nx_graph, _network.__getattribute__(f"APST_{4}") ,
                                        Cap_core, Cap_access, V, Nexullance_MP_APST4_M_s[1])
    nexu.init_model()
    _, MP_APST4_RT = nexu.solve()
    for M_EPs in M_EPs_s:
        # apply this routing table on all traffic demand matrices
        core_link_flows, access_link_flows = _network.distribute_M_EPs_on_weighted_paths(MP_APST4_RT, EPR, M_EPs)
        # calculate phi
        Nexullance_MP_APST4_shift_1_phis.append(gl.network_total_throughput(M_EPs, max(core_link_flows)/Cap_core, max(access_link_flows)/Cap_access)/(V*EPR))
    # =================================================================

    # calculate the routing table of Nexullance_MP_APST4 for 'half_shift'
    nexu = Nexullance_MP.Nexullance_MP(_network.nx_graph, _network.__getattribute__(f"APST_{4}") ,
                                        Cap_core, Cap_access, V, Nexullance_MP_APST4_M_s[2])
    nexu.init_model()
    _, MP_APST4_RT = nexu.solve()
    for M_EPs in M_EPs_s:
        # apply this routing table on all traffic demand matrices
        core_link_flows, access_link_flows = _network.distribute_M_EPs_on_weighted_paths(MP_APST4_RT, EPR, M_EPs)
        # calculate phi
        Nexullance_MP_APST4_halfshift_phis.append(gl.network_total_throughput(M_EPs, max(core_link_flows)/Cap_core, max(access_link_flows)/Cap_access)/(V*EPR))
    # =================================================================

    # calculate the phis with ECMP_ASP:
    for M_EPs in M_EPs_s:
        # apply this routing table on all traffic demand matrices
        core_link_flows, access_link_flows = _network.distribute_M_EPs_on_weighted_paths(_network.ECMP_ASP, EPR, M_EPs)
        # calculate phi
        ECMP_ASP_phis.append(gl.network_total_throughput(M_EPs, max(core_link_flows)/Cap_core, max(access_link_flows)/Cap_access)/(V*EPR))
    # =================================================================

    # now start with MD_Nexullance:
    MD_container = nexullance_exp_container(topo_name, V, D, EPR)
    MD_MP_obj, Nexullance_MD_MP_APST4_phis= MD_container.run_MD_nexullance_MP(M_EPs_s, M_weights, 4)
    MD_IT_obj, Nexullance_MD_IT_phis= MD_container.run_MD_nexullance_IT(M_EPs_s, M_weights)

    # now apply Nexullance_MP_APST4 on each individual traffic demand matrices:
    for i, M_EPs in enumerate(M_EPs_s):
        MP_phi, _ = MD_container.run_nexullance_MP(4, M_EPs, M_EPs_names[i])
        Nexullance_MP_APST4_phis.append(MP_phi)
        MP_phi, _ = MD_container.run_nexullance_IT(M_EPs, M_EPs_names[i])
        Nexullance_IT_phis.append(MP_phi)
    # =================================================================


    # calculate the objective function of each method:
    print("Objective function := the harmonic mean of network data throughout under demand matrices")
    # ECMP_ASP:
    print("ECMP_ASP objective function:", gl.cal_MD_obj_func(ECMP_ASP_phis, M_weights))
    # Nexullance_MP_APST4:
    print("Nexullance_MP_APST4 objective function:", gl.cal_MD_obj_func(Nexullance_MP_APST4_phis, M_weights))
    # Nexullance_IT:
    print("Nexullance_IT objective function:", gl.cal_MD_obj_func(Nexullance_IT_phis, M_weights))
    # MD_Nexullance_MP_APST4:
    assert(round(MD_MP_obj, 2) == round(gl.cal_MD_obj_func(Nexullance_MD_MP_APST4_phis, M_weights), 2))
    print("MD_Nexullance_MP_APST4 objective function:", gl.cal_MD_obj_func(Nexullance_MD_MP_APST4_phis, M_weights))
    # MD_Nexullance_IT:
    assert(round(MD_IT_obj, 2) == round(gl.cal_MD_obj_func(Nexullance_MD_IT_phis, M_weights), 2))
    print("MD_Nexullance_IT objective function:", gl.cal_MD_obj_func(Nexullance_MD_IT_phis, M_weights))
    # Nexullance_MP_APST4_uniform
    print("Nexullance_MP_APST4_uniform objective function:", gl.cal_MD_obj_func(Nexullance_MP_APST4_uniform_phis, M_weights))
    # Nexullance_MP_APST4_shift_1
    print("Nexullance_MP_APST4_shift_1 objective function:", gl.cal_MD_obj_func(Nexullance_MP_APST4_shift_1_phis, M_weights))
    # Nexullance_MP_APST4_halfshift
    print("Nexullance_MP_APST4_halfshift objective function:", gl.cal_MD_obj_func(Nexullance_MP_APST4_halfshift_phis, M_weights))


    filename = f'{topo_name}_{V}_{D}_MD_IT_uniform_shift_1_halfshift.csv'
    # save data to csv file
    with open(filename, 'a', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(csv_cols)
        csvfile.flush()
        for i in range(len(M_EPs_s)):
            csvwriter.writerow([V, D, M_EPs_names[i], M_weights[i], ECMP_ASP_phis[i], Nexullance_MP_APST4_phis[i],
                                Nexullance_IT_phis[i], Nexullance_MD_IT_phis[i], Nexullance_MD_MP_APST4_phis[i],
                                Nexullance_MP_APST4_uniform_phis[i], Nexullance_MP_APST4_shift_1_phis[i],
                                Nexullance_MP_APST4_halfshift_phis[i]
                                ])
            csvfile.flush()


if __name__ == '__main__':
    main()