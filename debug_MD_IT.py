import sys
sys.path.append("/users/ziyzhang/topology-research/nexullance/IT_boost/build")
from Nexullance_IT_cpp import MD_Nexullance_IT_interface

sys.path.append("/users/ziyzhang/topology-research")
import globals as gl
import topologies.RRG as RRG
import numpy as np

V = 16
D = 5
EPR = (D+1)//2
_network = RRG.RRGtopo(V, D)
ASP, _ = _network.calculate_all_shortest_paths()
ECMP_ASP = gl.ECMP(ASP)
arcs = _network.generate_graph_arcs()

Cap_remote = 10 #GBps
Cap_local = 10 #GBps

M_EPs = gl.generate_uniform_traffic_pattern(V, EPR)
remote_link_flows, local_link_flows = _network.distribute_M_EPs_on_weighted_paths(ECMP_ASP, EPR, M_EPs)
max_remote_link_load = np.max(remote_link_flows)/Cap_remote
max_local_link_load = np.max(local_link_flows)/Cap_local
# adapt the traffic scaling factor to 10x saturation
traffic_scaling = 10.0/max(max_local_link_load, max_remote_link_load)
M_EPs = traffic_scaling * M_EPs
remote_link_flows, local_link_flows = _network.distribute_M_EPs_on_weighted_paths(ECMP_ASP, EPR, M_EPs)
M_R = gl.convert_M_EPs_to_M_R(M_EPs, V, EPR)

nexu_it = MD_Nexullance_IT_interface(V, arcs, [M_R], [0.5], False)
# nexu_it.set_parameters(0.1, 7.0)
nexu_it.run()
nexu_it.get_max_link_load()