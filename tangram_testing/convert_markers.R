#  This R script is paired with the analysis done by the scripts
#  "four_sample_demo.py" and "four_sample_demo.sh". It reformats the marker
#  gene list produced by the script local_tangram_expdata/sce_to_anndata.R.

library('here')

orig_marker_path = file.path(here::here(), 'local_tangram_expdata', 'out', 'markers.csv')
new_marker_path = file.path(here::here(), 'tangram_testing', 'markers.txt')

#  These genes were named individually by Kristen over slack, and should be
#  used to test the Tangram mapping (not train on!)
test_markers = c('SNAP25', 'MBP', 'PCP4', 'CCK', 'RORB', 'ENC1', 'CARTPT', 
                 'NR4A2', 'RELN')
 
#  Read in any clean up marker list
markers = read.table(orig_marker_path)
markers = gsub('"|,', '', markers[2:nrow(markers), 2])

#  Confirm that none of the test markers named by Kristen are in the list of
#  genes to train on
if (any(test_markers %in% toupper(markers)) {
    stop('None of the test markers should be in the training genes!')
}

writeLines(markers, con=new_marker_path)
