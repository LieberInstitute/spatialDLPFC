#   Originally, SPOTlight was run in a way that was not reproducible, and the
#   SCE object was susbet to only 100 cells per cell type, which was suspected
#   to lead to poor results. This script keeps a copy of the old results (for
#   comparison) and moves the newer reproducible results to the original
#   location for use downstream. This was run interactively.

base_dir=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
good_results_name=subset_n1000

for dataset in IF nonIF; do
    for resolution in broad layer; do
        echo $dataset $resolution
        plot_dir=$base_dir/plots/spot_deconvo/04-spotlight/$dataset/$resolution

        #   Make a copy of the old plots
        mkdir $plot_dir/subset_n100_bad
        cp $plot_dir/*.pdf $plot_dir/subset_n100_bad/

        #   Move the good results to the old location
        cp $plot_dir/$good_results_name/*.pdf $plot_dir/ 
    done
done
