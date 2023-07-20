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
* To relate original file names to files shown here use _TODO by Nick_. Example: _TODO by Nick_.

### snRNA-seq

FASTQs from single-nucleus sequencing (those in `snRNA-seq`) use the following naming convention:

`[donor]_[position]_[mate]_[file number]` (e.g. `Br2720_mid_R1_1.fastq.gz`)

As with the spatial FASTQs, `[file number]` is added arbitrarily to distinguish separate files used to sequence the same sample and mate, and `[donor]` and `[position]` reflect the true sample ID (after identity swaps).

Related files:

* [`spaceranger` script](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/spaceranger_IF/01_spaceranger_IF.sh) + [`spaceranger` input parameters](https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/spaceranger_IF/spaceranger_IF_parameters.txt) (using original file names).
* Sample information [file](https://github.com/LieberInstitute/spatialDLPFC/blob/main/raw-data/sample_info/Visium_IF_DLPFC_MasterExcel_01262022.xlsx)
* To relate original file names to files shown here use _TODO by Nick_. Example: _TODO by Nick_.

## Relevant code

See `code/synapse_upload/01-prepare_fastq/01-prepare_fastq.R` and `code/globus/01-fastq_copy.R` (run in that order) for exact code used to copy and rename FASTQs from their original locations and naming conventions. Note that in the former script, `code/spaceranger/spaceranger_parameters.txt` was updated to use FASTQ filenames with the new naming convention (effectively resolving sample-ID mismatches as described above).

Related files:

* [`cellranger` scripts](https://github.com/LieberInstitute/DLPFC_snRNAseq/tree/main/code/01_align) (using original file names).
* Sample information [files](https://github.com/LieberInstitute/DLPFC_snRNAseq/tree/main/raw-data/sample_info) extracted to a [smaller csv file](https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/processed-data/04_synapse_upload/pd.csv)
* [Table](https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/processed-data/04_synapse_upload/fastq_renaming_scheme.csv) showing how original file names related to renamed file names uploaded to Synapse (same files shown here in `FASTQ_globus`).
* To relate original file names to files shown here use _TODO by Nick_. Example: _TODO by Nick_.
