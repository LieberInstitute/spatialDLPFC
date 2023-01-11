#   Originally, SPOTlight was run in a way that was not reproducible, and the
#   SCE object was susbet to only 100 cells per cell type, which was suspected
#   to lead to poor results. This script keeps a copy of the old results (for
#   comparison) and moves the newer reproducible results to the original
#   location for use downstream. This was run interactively.

base_dir=/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC
good_results_name=full_data

for dataset in IF nonIF; do
    for resolution in broad layer; do
        plot_dir=$base_dir/plots/spot_deconvo/04-spotlight/$dataset/$resolution
        processed_dir=$base_dir/processed-data/spot_deconvo/04-spotlight/$dataset/$resolution

        ########################################################################
        #   Make a copy of the old results
        ########################################################################

        #   Create locations to copy the old results to
        mkdir $plot_dir/subset_n100_bad
        mkdir $processed_dir/subset_n100_bad

        #   Keep just the PDF plots
        cp $plot_dir/*.pdf $plot_dir/subset_n100_bad/

        #   Copy the sample folders containing 'clusters.csv' files and the RDS
        #   results (everything except the 'subset_n*' and 'full_data' folders)
        cp -R $processed_dir/Br*_* $processed_dir/subset_n100_bad/
        cp $processed_dir/*.rds $processed_dir/subset_n100_bad/

        ########################################################################
        #   Move the good results to the old location
        ########################################################################
        
        cp $plot_dir/$good_results_name/*.pdf $plot_dir/
        cp -R $processed_dir/$good_results_name/* $processed_dir/
    done
done
