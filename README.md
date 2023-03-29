# Diagnotics and Analysis Codes for MIDAS2

## igv_render (2023-03-11)

Modified based on [sbx_igv](https://github.com/sunbeam-labs/sbx_igv), igv_render generate an alignment image, given a genome, list of bam files, list of segments. IGVtools commands can refer to [here](https://software.broadinstitute.org/software/igv/automation).

### Install IGVtools and pysam

```
conda install -c bioconda igv
conda install -c conda-forge xvfbwrapper
conda install -c bioconda pysam
pip install python-libxdo
```

### Example Command
```
python render_igv.py \
    --genome_fasta /pollard/home/czhao/igv_data/bt2_index_biohub/repgenomes.fa \
    --bamfiles example_data/list_of_bams  --segment_file example_data/list_of_segments \
    --outdir example_data/igv_outputs 
```



