# Modified from https://github.com/sunbeam-labs/sbx_igv
# Chunyu Zhao 2023-03-10

"""
Helper functions for interfacing with the Integrative Genomics Viewer (IGV).

Dependencies:
 * IGV
 * xvfb
 * xdotool (for the socket-based method)
"""

import socket
import subprocess
import tempfile
import time
import os
import argparse
import pysam
from pathlib import Path
from Bio import SeqIO
from utils import command
from collections import defaultdict


IGV_PREFS = {"ENABLE_ANTIALIASING":True, "NAME_PANEL_WIDTH":200, "SAM.MAX_VISIBLE_RANGE": 1000}


def render(genome, bams, imagefile, seqID=None, igv_fp="igv", method="script", locus = None, igv_prefs=None):
        """ Render an alignment to an image, given a genome and bam files.
        genome: path to a fasta file
        bams: list of path to a sorted, indexed bam file
        imagefile: path to the image to save
        seqID: (optional) sequence identifier to load from genome
        igv_fp: (optional) path to IGV executable
        method: (optional) method for controlling IGV process: "script" or "socket"
        igv_prefs: (optional) dictionary of IGV preferences to override

        The image file may be smaller than expected.  See
        igv_render_socket_nonblocking() for an attempt to enlarge the window
        before saving the image.
        """
        igv_prefs = igv_prefs or {}
        genome = str(genome)
        imagefile = str(imagefile)
        input_paths = [str(Path(str(bam)).resolve()) for bam in bams]
        genome_path = str(Path(genome).resolve())
        output_path = str( Path('.').resolve() / Path(imagefile) )
        # build a "seqID:1:length" string to force IGV to display the full
        # segment.  If no segment was given it will default to the first.
        if len(locus):
            goto_locus = locus
        else:
            goto_locus = _seq_length(genome, seqID)
            goto_locus = "%s:%s-%s" % (goto_locus[0], '1', goto_locus[1])
        igvcommands = ['new',
            'genome ' + genome_path,
            'goto ' + goto_locus,
            'load ' + ','.join(input_paths),
            'collapse',
            'snapshot ' + output_path,
            'exit']
        # If previous genome files listed in IGV's preferences are no
        # longer available, IGV will throw a null pointer exception at
        # startup and batch commands will fail.  So, we'll use a
        # preferences override file to list the genome file used here.
        # (This is probably an IGV bug.  We should see if it happens in the
        # latest release.) I've also tried setting IGV.Bounds in an attempt
        # to make the window larger, but it doesn't seem to have any
        # effect.
        igv_prefs["GENOME_LIST"] = ";" + genome_path
        igv_prefs["DEFAULT_GENOME_KEY"] = genome_path
        print(igvcommands)
        _control_script(igvcommands, igv_fp, igv_prefs)


def _control_script(igvcommands, igv_fp, igv_prefs):
        igvscript = tempfile.NamedTemporaryFile()
        igvscript.writelines(map(lambda x: bytes(x+'\n', 'ascii'), igvcommands))
        igvscript.flush()
        igvprefsfile = _write_prefs(igv_prefs)
        command("xvfb-run -a -s '-screen 1 1920x1080x24' %s -o %s -b %s" % (igv_fp, igvprefsfile.name, igvscript.name))


def _write_prefs(igv_prefs):
        igvprefsfile = tempfile.NamedTemporaryFile()
        for k,v in igv_prefs.items():
            igvprefsfile.write(bytes("%s=%s\n" % (k, v), 'ascii'))
        igvprefsfile.flush()
        return igvprefsfile


def _seq_length(fasta_fp, seqID=None):
    """Give the ID and sequence length of the first (or specified) sequence. """
    record_lengths = [[r.id, len(r.seq)] for r in SeqIO.parse(fasta_fp, "fasta")]
    if not seqID:
        return(record_lengths[0])
    record_length_match = [r for r in record_lengths if r[0] == seqID]
    return(record_length_match[0])


def read_seq_ids(fasta_fp):
    """
    Return the sequence identifiers for a given fasta filename.
    """
    return [record.id for record in SeqIO.parse(str(fasta_fp), "fasta")]


def read_bamlist(bamfile_fp):
    list_of_bams = []
    with open(bamfile_fp) as stream:
        for line in stream:
            list_of_bams.append(line.strip("\n"))
    return list_of_bams


def load_locus_list(gseg_fp):
    list_of_locus = []
    with open(gseg_fp) as stream:
        for line in stream:
            list_of_locus.append(line.strip("\n"))
    return list_of_locus


def compute_genome_length(genome_file):
    genome_len_dict = defaultdict(int)
    with open(genome_file) as stream:
        for sn, rec in enumerate(SeqIO.parse(stream, 'fasta')):
            contig_seq = str(rec.seq).upper()
            contig_len = len(contig_seq)
            genome_len_dict[rec.id] = contig_len
    return genome_len_dict


def main():
    p = argparse.ArgumentParser(prog="python summarize_pileup.py", description='compute sites level summary per pileup file')
    p.add_argument(
        "--genome_fasta", required=True,
        type=str,
        help=f"Path of MIDAS2 repgenome fasta file.")
    p.add_argument(
        "--bamfiles", required=True,
        type=str,
        help=f"List of path to sorted and indexed BAM files.")
    p.add_argument(
        "--outdir", required=True,
        type=str,
        help=f"Path to the base directory of rendered images.")
    p.add_argument(
        "--segment_file", required=False,
        type=str,
        default = [],
        help=f"Genome segments. Format: CHROM:START-END. Per segment each row. CHROM is required.")


    args = p.parse_args()

    genome_file = args.genome_fasta
    glen_dict = compute_genome_length(genome_file)

    if len(args.segment_file):
        genome_segments = load_locus_list(args.segment_file)
        assert len(genome_segments), f"Genome Segment File Is Empty"
    else:
        print("This takes a long time")
        genome_segments = read_seq_ids(genome_file)

    list_of_bams = read_bamlist(args.bamfiles)

    Path(args.outdir).mkdir(parents = True, exist_ok=True)

    if ":" in genome_segments[0]:
        for locus in genome_segments:
            segment = locus.split(":")[0]
            outpng = f"{args.outdir}/{locus}.png"
            render(genome=genome_file, bams=sorted(list_of_bams), imagefile=outpng, seqID=segment, locus = locus, igv_prefs=IGV_PREFS)
    else:
        window_size = 40000
        for segment in genome_segments:
            segment_length = glen_dict[segment]
            for lower_bound in range(0, segment_length, window_size):
                upper_bound = min(segment_length, lower_bound + window_size)
                locus = f"{segment}:{lower_bound}-{upper_bound}"
                outpng = f"{args.outdir}/{locus}.png"
                render(genome=genome_file, bams=sorted(list_of_bams), imagefile=outpng, seqID=segment, locus = locus, igv_prefs=IGV_PREFS)

main()
