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

def main():
    config = gl.ddf_configs[0]
    V = config[0]
    D = config[1]
    EPR = (D+1)//2
    topo_name="DDF"
    _network = HPC_topo.HPC_topo.initialize_child_instance(topo_name+"topo", V, D)

    _network.pre_calculate_ECMP_ASP()

    # ==============================================================
    M_EPs_names=['uniform', 'shift_1', 'half_shift'] # traffic demand matrices to test on
    M_EPs_s = gl.gen_M_EPs_s(topo_name, V, D, EPR, M_EPs_names, 10.0)
    M_weights = [1/3, 1/3, 1/3]
    # M_weights = [0.1, 0.6, 0.3]

    # now start with MD_Nexullance:
    MD_container = nexullance_exp_container(topo_name, V, D, EPR)
    _result = MD_container.run_and_profile_MD_nexullance_IT(M_EPs_s, M_weights, 10)

    print(f"ave time for MD_IT is {_result["ave_time[s]"]} s, std = {_result["std_time[s]"]} s")

    # MD_IT_obj, Nexullance_MD_IT_phis= MD_container.run_MD_nexullance_IT(M_EPs_s, M_weights, False)

    # for i, M_EPs in enumerate(M_EPs_names):
    #     print(f"phi for {M_EPs} is {Nexullance_MD_IT_phis[i]}")

    # # MD_Nexullance_IT:
    # assert(round(MD_IT_obj, 2) == round(gl.cal_MD_obj_func(Nexullance_MD_IT_phis, M_weights), 2))
    # print("MD_Nexullance_IT objective function:", gl.cal_MD_obj_func(Nexullance_MD_IT_phis, M_weights))


if __name__ == '__main__':
    main()