#$ -cwd
#$ -o ../../processed-data/04_profiling/01-squidpy_segment_single_core.log
#$ -e ../../processed-data/04_profiling/01-squidpy_segment_single_core.log
#$ -N squidpy_segment_single_core
#$ -l caracol,mf=10G,h_vmem=10G

module load tangram/1.0.2

echo "Detected $(nproc) cores usable of $(nproc --all) total cores on the node."

thread_list=(1 2 4 8 16)
for nthreads in "${thread_list[@]}"; do
    echo "Setting all environment variables to $nthreads..."
    export OMP_NUM_THREADS=$nthreads
    export MKL_NUM_THREADS=$nthreads
    export OPENBLAS_NUM_THREADS=$nthreads

    python 01-squidpy_segment.py
done
