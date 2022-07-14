library("styler")
library("biocthis")
library("here")

#  Style R scripts consistently. This script is intended to
#  be invoked whenever one or more R scripts is modified during development
style_dir(
    path = here("tangram_libd", "code", "03_nn_run"),
    transformers = bioc_style(),
    recursive = FALSE
)
