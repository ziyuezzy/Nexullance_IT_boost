cd FlameGraph
sudo perf record -g -F 100000 bash ../run.sh 
sudo perf script | ./stackcollapse-perf.pl > out.perf-folded
./flamegraph.pl out.perf-folded > ../flamegraph.svg