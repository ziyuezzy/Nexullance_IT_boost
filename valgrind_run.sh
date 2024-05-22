export LD_LIBRARY_PATH=/users/ziyzhang/topology-research/nexullance/IT_boost/boost/stage/lib

# ./main /users/ziyzhang/topology-research/nexullance/handover_data/RRG_16_5/graph.graphml
topo=RRG_16_5
# topo=RRG_100_11
traffic=uniform
valgrind --leak-check=full --show-leak-kinds=all /users/ziyzhang/topology-research/nexullance/IT_boost/main.exe /users/ziyzhang/topology-research/nexullance/handover_data/"$topo"/graph.graphml /users/ziyzhang/topology-research/nexullance/handover_data/"$topo"/"$traffic".txt 0
