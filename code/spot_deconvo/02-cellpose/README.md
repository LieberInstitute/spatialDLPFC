This directory contains work to use `cellpose` to segment nuclei on our IF data, ultimately with the goal of producing cell counts per spot (used for spot deconvolution, such as as a required input to `tangram`, and beneficial one to `cell2location`). An additional goal sought with the code here is take DAPI segmentations produced by `cellpose` and directly use fluorescence images to call cell types.

Old/ deprecated code:

* `02-count_cells.*`: directly adapting [Richard's code](https://github.com/chaichontat/libd-rotation/blob/main/scripts/segmentation/process_mask.py) to count cell types
* `03-dilate_masks.*`: experiment with a variation on image dilation to enlarge each nucleus mask, since we were interested in quantifying fluorescence in the nucleus as well as a small region around it
* `05-classify_nuclei.*`: a manual method simply taking the maximum (normalized) fluoresence in each image channel to call cell types. This method was more complicated and performed worse than it's replacement: a `DecisionTreeClassifier`, used in the scripts `07-cart.*` and `08_classify_nuclei_cart.*`
* `06-evaluate_method.py`: used to assess the accuracy of the manual method used in `05-classify_nuclei.*`

Latest workflow (in order):

* `01-create_masks.*`: segment nuclei on each IF image using `cellpose` to output masks
* `04-quantify_fluor.*`: quantify fluorescence within each nucleus (as well as a small area around it) to output a `DataFrame` of intensities for each cell
* `09-prepare_loopy.*`: prepare a CSV of cell IDs and coordinates supported as input to the [Loopy Browser](https://loopybrowser.com/) so that cells can be manually labelled by cell type (to train a model).
* `07-cart.*`: train a `DecisionTreeClassifier` on manually annotated cells, and save the model
* `08-classify_nuclei_cart.*`: Use the saved `DecisionTreeClassifier` to classify cell types for all IF samples and all cells. Output a CSV of metrics (rows are cells and columns are things like fluorescence intensity, cell area, etc) and CSV of cell counts (rows are spots and columns are cell-type counts)
