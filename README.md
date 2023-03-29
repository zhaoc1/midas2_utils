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

## snv_track (2023-03-09)

R script to generate SNV track plot based on MIDAS2 `snps_info.tsv`.

TODO: 
- Either keep the fixed sites for the `merge_snps`, Or compute the mean site depth across samples for the fixed sites.


## Blast results plot

```
df <- read_delim("blast_outfmt6.tsv", delim = "\t", col_names = F, show_col_types = F)

df <- df %>% set_colnames(c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "e_value", "bit_score"))

df %>%
  mutate(xmin = pmin(qstart, qend), xmax = pmax(qstart, qend)) %>%
  ggplot() + 
  geom_rect(aes(xmin = xmin, xmax = xmax, fill=pident), ymin = -Inf, ymax = Inf, color = NA, alpha = 0.9) + 
  #geom_text(aes(x = xmin, y = 0.01, label=baiGene), size = 3, vjust = 0, hjust = 0, check_overlap = FALSE) +
  geom_vline(aes(xintercept = as.numeric(12000)), colour = "red", alpha = 0.8) +
  #ylim(c(0, 0.1)) +
  theme_bw() +
  scale_fill_viridis(alpha=0.9, discrete=FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~sseqid, ncol = 1)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), xis.text.y=element_blank(), axis.ticks.y = element_blank())
```
