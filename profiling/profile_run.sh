for num_samples in 1 2 4 ; do
    sudo bash gen_flamegraph.sh profile_diff_IT.py --num-samples "$num_samples" --output FFT3D_"$num_samples"_sample_diff_IT --freq 99
    sudo bash gen_flamegraph.sh profile_IT.py --num-samples "$num_samples" --output FFT3D_"$num_samples"_sample_IT --freq 99
done