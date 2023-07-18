Copies of the raw FASTQ files for the entire project, accessible through our Globus endpoint.

Note that we also have the `raw-data/FASTQ` directory, which consists of symlinks and therefore more explicitly links back to the original FASTQs as uploaded after sequencing. The copies present in this directory are necessary to isolate this project's FASTQs to one parent directory, for ease of access through Globus (as opposed to symlinks that point to locations we aren't sharing).
