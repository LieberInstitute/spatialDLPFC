library("here")

## Check sample_ids from
## https://github.com/LieberInstitute/spatialDLPFC/blob/4c127dafb03fa1b185168a22544453b55b833d97/code/analysis/01_build_SPE_final.R#L35-L67
sample_id <- c(
    "DLPFC_Br2743_ant_manual_alignment",
    "DLPFC_Br2743_mid_manual_alignment_extra_reads",
    "DLPFC_Br2743_post_manual_alignment",
    "DLPFC_Br3942_ant_manual_alignment",
    "DLPFC_Br3942_mid_manual_alignment",
    "DLPFC_Br3942_post_manual_alignment",
    "DLPFC_Br6423_ant_manual_alignment_extra_reads",
    "DLPFC_Br6423_mid_manual_alignment",
    "DLPFC_Br6423_post_extra_reads",
    "DLPFC_Br8492_ant_manual_alignment",
    "DLPFC_Br8492_mid_manual_alignment_extra_reads",
    "DLPFC_Br8492_post_manual_alignment",
    "Round4/DLPFC_Br2720_ant_2",
    "Round2/DLPFC_Br2720_mid_manual_alignment",
    "Round2/DLPFC_Br2720_post_extra_reads",
    "Round4/DLPFC_Br6432_ant_2",
    "Round2/DLPFC_Br6432_mid_manual_alignment",
    "Round2/DLPFC_Br6432_post_manual_alignment",
    "Round3/DLPFC_Br6471_ant_manual_alignment_all",
    "Round3/DLPFC_Br6471_mid_manual_alignment_all",
    "Round3/DLPFC_Br6471_post_manual_alignment_all",
    "Round3/DLPFC_Br6522_ant_manual_alignment_all",
    "Round3/DLPFC_Br6522_mid_manual_alignment_all",
    "Round3/DLPFC_Br6522_post_manual_alignment_all",
    "Round3/DLPFC_Br8325_ant_manual_alignment_all",
    "Round4/DLPFC_Br8325_mid_2",
    "Round3/DLPFC_Br8325_post_manual_alignment_all",
    "Round3/DLPFC_Br8667_ant_extra_reads",
    "Round3/DLPFC_Br8667_mid_manual_alignment_all",
    "Round3/DLPFC_Br8667_post_manual_alignment_all"
)

df <- read.table(here("code", "spaceranger", "spaceranger_parameters.txt"))
stopifnot(all(basename(sample_id) %in% df$V1))
