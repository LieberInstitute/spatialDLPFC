This directory contains work done to prepare directories compatible with the `Loopy Browser`, where we'll ultimately host the 4 IF and 30 H&E images, with spot deconvolution results as displayable features in the browser.

- `01-host_image_local.sh`: first test, creating the files for Loopy to host just the IF image for one sample.
- `02-IF.*`: a more complete re-writing of `01-host_image_local.sh` that produces directories for locally hosting all 4 IF samples and their corresponding spot-deconvolution results. Specifically, software-estimated results are provided at broad and layer resolution, and CART results are provided at collapsed resolution.
- `03-nonIF.*`: a corresponding version of `02-IF.*` for the n=30 nonIF (H&E) samples. Software-estimated results are provided at broad and layer resolution.
