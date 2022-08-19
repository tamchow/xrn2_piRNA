import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

from get_seq import get_seq

samples = sys.argv[1].split(":")
query_pirna = sys.argv[2]
prefix = sys.argv[3]
titles = sys.argv[4].split(":") if len(sys.argv) > 4 else samples

min_overlap_frac = 0.3

transcripts = pd.read_csv(
    "ref/mirnas_pirnas.tsv",
    sep="\t",
    index_col="transcript_id",
    dtype={
        "transcript_id": str,
        "Chr": "category",
        "Start": int,
        "End": int,
        "Strand": "category",
        "Length": int,
        "gene_id": str,
        "biotype": "category",
        "name": str,
    },
)

transcripts = transcripts[~transcripts.index.duplicated(keep="first")]

transcripts["Start"] = transcripts["Start"] - 1
transcripts["End"] = transcripts["End"]

if "name" not in transcripts.columns:
    transcripts["name"] = transcripts.index.copy()


def reads_mapped_to_region(alns, chromosome, start, end):
    region_len = abs(end - start)
    min_overlap = min_overlap_frac * region_len
    aln_count = 0
    for aln in alns.fetch(chromosome, start, end):
        overlap = aln.get_overlap(start, end)
        if overlap is not None and overlap >= min_overlap:
            aln_count += 1
    return aln_count


def coverage_hist(pileup, start, end):
    return np.array([pileup_col.get_num_aligned() for pileup_col in pileup if start <= pileup_col.reference_pos <= end])


def coverage(sample, chromosome, start, end, scale_factor=1):
    alns = pysam.AlignmentFile(f"{prefix}/{sample}.sorted.bam", "rb")
    raw_hist = coverage_hist(alns.pileup(chromosome, start, end + 1, truncate=False), start, end)
    mapped_reads = reads_mapped_to_region(alns, chromosome, start, end)
    return raw_hist / mapped_reads * scale_factor


upstream_nt, feature_nt, total_nt = 2, 21, 30

fig = plt.figure(figsize=(12, 3 * len(samples)), dpi=300)

axs = fig.subplots(nrows=len(samples), ncols=1)
for ax, sample, title in zip(axs, samples, titles):
    try:
        query_row = transcripts[transcripts['name'] == query_pirna]
        chrom, start, end, strand = (
            query_row["Chr"].values[0],
            query_row["Start"].values[0],
            query_row["End"].values[0],
            query_row["Strand"].values[0]
        )
        print(query_pirna, chrom, start, end, strand)
        downstream_nt = total_nt - (upstream_nt + feature_nt)
        sequence = get_seq('ref/c_elegans.PRJNA13758.WS280.genomic.fa', chrom, start + 1, end, strand, upstream_adjust=upstream_nt, downstream_adjust=downstream_nt)  # get_seq expects 1-based coordinates
        start = start - upstream_nt
        end = start + total_nt
        print(query_pirna, chrom, start, end, strand, sequence)
        x = list(range(-upstream_nt, 0)) + list(range(1, total_nt - upstream_nt + 1))  # since there is no 0th coordinate
        x = [f"{idx}\n{base}" for idx, base in zip(x, sequence)]
        y = coverage(sample, chrom, start, end)
        y = np.concatenate((y[:len(x)], np.zeros(max(len(x) - len(y), 0))))  # ensure y is the same length as x
        tmp_df = pd.DataFrame(data={'Base Position': x, 'Normalized Coverage': y})
        sns.barplot(data=tmp_df, x='Base Position', y='Normalized Coverage', ax=ax, color='tab:blue')
        ax.set_ylim(ymin=0, ymax=1)
        ax.set_title(title)
    except (KeyError, IndexError):
        print(f'Could not find {query_pirna} in the list of known piRNAs - please check the format of 21ur-xxxx or piR-cel-xxxx')

fig.subplots_adjust(hspace=0.8)
fig.savefig(f"{prefix}/{query_pirna}_coverage.png")
