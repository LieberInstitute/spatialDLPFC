Copies of the raw FASTQ files for the entire project, accessible through our Globus endpoint.

## Motivation

Note that we also have the `raw-data/FASTQ` directory, which consists of symlinks and therefore more explicitly links back to the original FASTQs as uploaded after sequencing. The copies present in this directory are necessary to isolate this project's FASTQs to one parent directory, for ease of access through Globus (as opposed to symlinks that point to locations we aren't sharing). This is due to restrictions Globus has on symlinks described in [their documentation](https://docs.globus.org/faq/transfer-sharing/#how_does_globus_handle_symlinks).

## Naming convention

[Renaming script](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/globus/01-fastq_copy.R).

### Visium/ Visium-SPG

The spatial FASTQs (those in `Visium`/`Visium-SPG`) use the following naming convention:

`[slide]_[array]_[mate]_[file number].fastq.gz` (e.g. `V10U24-091_A1_R1_1.fastq.gz`)

Here `[file number]` is necessary because the same sample and mate may be split across multiple lanes. Note that this naming convention discards sequencing lane, and `[file number]` is chosen arbitrarily simply to make filenames unique.

Note also that in some cases, samples IDs were swapped after discovering mislabelled samples. The file names here reflect the true sample IDs, and thus some FASTQs here contain data originally matching FASTQs named with a different sample ID.

Related files:

* [`spaceranger` script](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/spaceranger/rerun_spaceranger.sh) + [`spaceranger` input parameters](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/spaceranger/spaceranger_parameters.txt) (using original file names).
* Sample information [file](https://github.com/LieberInstitute/spatialDLPFC/blob/main/raw-data/sample_info/Visium_dlpfc_mastersheet.xlsx)
* To relate original file names to files shown here use [this table](https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/raw-data/FASTQ_Globus/fastq_mapping.csv), described more thoroughly [below](https://github.com/LieberInstitute/spatialDLPFC/tree/main/raw-data/FASTQ_Globus#relevant-code)

### snRNA-seq

FASTQs from single-nucleus sequencing (those in `snRNA-seq`) use the following naming convention:

`[donor]_[position]_[mate]_[file number]` (e.g. `Br2720_mid_R1_1.fastq.gz`)

As with the spatial FASTQs, `[file number]` is added arbitrarily to distinguish separate files used to sequence the same sample and mate, and `[donor]` and `[position]` reflect the true sample ID (after identity swaps).

Related files:

* [`spaceranger` script](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/spaceranger_IF/01_spaceranger_IF.sh) + [`spaceranger` input parameters](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/spaceranger_IF/spaceranger_IF_parameters.txt) (using original file names).
* Sample information [file](https://github.com/LieberInstitute/spatialDLPFC/blob/main/raw-data/sample_info/Visium_IF_DLPFC_MasterExcel_01262022.xlsx)
* To relate original file names to files shown here use [this table](https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/raw-data/FASTQ_Globus/fastq_mapping.csv), described more thoroughly [below](https://github.com/LieberInstitute/spatialDLPFC/tree/main/raw-data/FASTQ_Globus#relevant-code)

## Relevant code

See `code/synapse_upload/01-prepare_fastq/01-prepare_fastq.R` and `code/globus/01-fastq_copy.R` (run in that order) for exact code used to copy and rename FASTQs from their original locations and naming conventions. Note that in the former script, `code/spaceranger/spaceranger_parameters.txt` was updated to use FASTQ filenames with the new naming convention (effectively resolving sample-ID mismatches as described above).

Related files:

* [`cellranger` scripts](https://github.com/LieberInstitute/DLPFC_snRNAseq/tree/main/code/01_align) (using original file names).
* Sample information [files](https://github.com/LieberInstitute/DLPFC_snRNAseq/tree/main/raw-data/sample_info) extracted to a [smaller csv file](https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/processed-data/04_synapse_upload/pd.csv)
* [Table](https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/raw-data/FASTQ_Globus/fastq_mapping.csv) mapping between FASTQs with the three naming conventions used at different points in the project: originals from sequencing, those uploaded to Synapse, and those shared through Globus. It contains the following columns:
    * `sample_id`: the sample identifier used for Globus FASTQs, as described in the [naming convention sections](https://github.com/LieberInstitute/spatialDLPFC/tree/main/raw-data/FASTQ_Globus#naming-convention) above
    * `donor`: the donor associated with this sample
    * `fastq_synapse`: FASTQ filepath uploaded to Synapse
    * `fastq_globus`: FASTQ filepath in this directory, shared through Globus
    * `fastq_sequencing`: FASTQ filepath output from sequencer (prior to any sample-identity swaps)
    * `assay`: Which type of experiment this FASTQ was associated with ("Visium", "Visium-SPG", or "snRNA-seq")
For example, if we have the FASTQ file `/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/raw-data/FASTQ_Globus/snRNA-seq/Br3942_ant_R1_1.fastq.gz` as shared through Globus, one can find this file in the `fastq_globus` column, and match it to the original sequencer file `/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-11-22_KMay110521/11c_k_L001_ds.b485101e8841432892fb04b1c1d9a3a6/11c-k_S8_L001_R1_001.fastq.gz` (found in the `fastq_sequencing` column). Note that paths in the `fastq_synapse` and `fastq_sequencing` columns are not shared through Globus, and exist simply to provide a record of filenames at different stages of the project.
