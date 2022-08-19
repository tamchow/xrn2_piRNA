import numpy as np
import pandas as pd
from scipy import stats
import pysam

import os
import sys

sample = sys.argv[1]
prefix = sys.argv[2] if len(sys.argv) > 2 else '.'
threads = int(os.environ.get("threads", len(os.sched_getaffinity(0)) // 2))
# mode should be one of 'pirna', 'mirna', 'other'
rna_mode = os.environ.get("rna_mode", 'pirna')

transcripts = None
if rna_mode == 'other':
    transcripts = pd.read_csv("ref/snoRNAs_rRNAs.tsv", sep='\t',
                              index_col='transcript_id',
                              dtype={'transcript_id': str,
                                     'Chr': 'category',
                                     'Start': int,
                                     'End': int,
                                     'Strand': 'category',
                                     'Length': int,
                                     'gene_id': str,
                                     'gene_biotype': 'category'
                                     })
else:
    transcripts = pd.read_csv("ref/mirnas_pirnas.tsv", sep='\t',
                              index_col='transcript_id',
                              dtype={'transcript_id': str,
                                     'Chr': 'category',
                                     'Start': int,
                                     'End': int,
                                     'Strand': 'category',
                                     'Length': int,
                                     'gene_id': str,
                                     'biotype': 'category',
                                     'name': str,
                                     })

transcripts = transcripts[~transcripts.index.duplicated(keep='first')]

transcripts['Start'] = transcripts['Start'] - 1  # start: convert from 1-based, inclusive to 0-based, inclusive
transcripts['End'] = transcripts['End']  # end: convert from 1-based, inclusive to 0-based, exclusive (i.e., leave it as-is)

if 'name' not in transcripts.columns:
    transcripts['name'] = transcripts.index.copy()

sam_file = pysam.AlignmentFile(os.path.join(prefix, f"{sample}.sorted.bam"), 'rb', threads=threads)

# miRDeep2 default
min_read_len, max_read_len = 18, 40
upstream_bases, downstream_bases = 2, 5

# # lenient
# min_read_len, max_read_len = 13, None
# upstream_bases, downstream_bases = 2, 5

# # stringent
# min_read_len, max_read_len = 20, None
# upstream_bases, downstream_bases = 2, 3

if rna_mode == 'other':
    # not small RNAs
    min_read_len, max_read_len = None, None
    upstream_bases, downstream_bases = 0, 0

check_strand = True

# Fraction of read length; 85% overlap to define piRNA loci or 99% to define miRNA loci, 10% overlap to define other loci
min_overlap_fracs = {'other': 0.10, 'pirna': 0.85, 'mirna': 0.99}
min_overlap_frac = min_overlap_fracs.get(rna_mode, min_overlap_fracs['pirna'])


def get_strand(read):
    if read.is_reverse:
        return '-'
    return '+'


def is_valid_read(read, start, end, strand):
    overlap = read.get_overlap(start, end)
    read_len = read.query_length
    min_overlap = min_overlap_frac * read.query_alignment_length
    if rna_mode == 'other':
        return False if overlap is None else overlap >= min_overlap
    return all((False if overlap is None else overlap >= min_overlap,
                True if min_read_len is None else read_len >= min_read_len,
                True if max_read_len is None else read_len <= max_read_len,
                get_strand(read) == strand if check_strand else True))
                # read.reference_start >= max(start - upstream_bases, 0),
                # False if read.reference_end is None else read.reference_end <= end + downstream_bases))


processed_reads = {}


def transcript_read_distribution(transcript_id):
    transcript_data = transcripts.loc[transcript_id]
    read_starts, read_ends = [], []
    try:
        # print(transcript_id, *transcript_data.to_list(), file=sys.stderr)
        start, end, strand = transcript_data['Start'], transcript_data['End'], transcript_data['Strand']
        for read in sam_file.fetch(transcript_data['Chr'], start, end):
            # print(f"{transcript_data['name']} -> {read.query_name}", file=sys.stderr)
            if not processed_reads.get(transcript_data['name']):
                processed_reads[transcript_data['name']] = set()
            if is_valid_read(read, start, end, strand):  # and read.query_name not in processed_reads[transcript_data['name']]:
                aln_start, aln_end = read.reference_start, read.reference_end
                # note aln_end is one past last aligned base
                read_starts.append(aln_start)
                read_ends.append(aln_end)

                processed_reads[transcript_data['name']].add(read.query_name)
                print(read.query_name, aln_start, aln_end, aln_end - aln_start, read.query_length)
    except ValueError:
        pass
    except Exception as e:
        print(f'Error {str(e)} processing transcript_id {transcript_id}', file=sys.stderr)
        raise e
    return tuple(read_starts), tuple(read_ends)


def safe_stdev(xs):
    if len(xs) > 1:
        return np.std(xs)
    return 0


def only_mode(xs):
    mode_, count_ = stats.mode(xs)
    return mode_.item()


def safe_statsfunc(func, xs):
    if len(xs) > 0:
        return func(xs)
    return np.nan


count_5p_ext, count_3p_ext = upstream_bases, downstream_bases
count_5p_ext_min_len, count_3p_ext_min_len = 23, 23


def read_pos_stats_(read_starts, ref_start, read_ends, ref_end):
    # 5p extension means mean_read_start < transcript_start
    read_starts, read_ends = np.array(read_starts), np.array(read_ends)
    # adjust end coordinates from 0-based, exclusive to 0-based, inclusive
    read_ends = read_ends - 1
    ref_end = ref_end - 1
    read_5p_exts, read_3p_exts = np.maximum(ref_start - read_starts, 0), np.maximum(read_ends - ref_end, 0)
    read_lengths = np.maximum(read_ends - read_starts, 0)
    start_stats = (safe_statsfunc(np.mean, read_5p_exts),
                   safe_statsfunc(safe_stdev, read_5p_exts),
                   safe_statsfunc(only_mode, read_5p_exts),
                   np.sum(read_5p_exts >= 1),
                   np.sum(read_5p_exts == count_5p_ext),
                   np.sum(read_5p_exts == 0),
                   np.sum(np.logical_and(read_5p_exts >= 1, read_lengths >= count_5p_ext_min_len)),
                   np.sum(np.logical_and(read_5p_exts == count_5p_ext, read_lengths >= count_5p_ext_min_len)),)
    end_stats = (safe_statsfunc(np.mean, read_3p_exts),
                 safe_statsfunc(safe_stdev, read_3p_exts),
                 safe_statsfunc(only_mode, read_3p_exts),
                 np.sum(read_3p_exts >= 1),
                 np.sum(read_3p_exts == count_3p_ext),
                 np.sum(read_3p_exts == 0),
                 np.sum(np.logical_and(read_3p_exts >= 1, read_lengths >= count_3p_ext_min_len)),
                 np.sum(np.logical_and(read_3p_exts == count_3p_ext, read_lengths >= count_3p_ext_min_len)),)
    len_stats = (safe_statsfunc(np.mean, read_lengths),
                 safe_statsfunc(safe_stdev, read_lengths),
                 safe_statsfunc(only_mode, read_lengths))
    return start_stats, end_stats, len_stats, np.sum(np.logical_and(read_5p_exts == 0, read_3p_exts == 0)), len(read_lengths)


def read_pos_stats(transcript_id):
    read_starts, read_ends = transcripts.loc[transcript_id]['read_distrib']
    return read_pos_stats_(read_starts, transcripts.loc[transcript_id]['Start'], read_ends, transcripts.loc[transcript_id]['End'])


transcripts['read_distrib'] = transcripts.index.map(transcript_read_distribution)
# transcripts.to_pickle(f'{prefix}/{sample}.pkl')
transcripts['5p_ext_stats'], transcripts['3p_ext_stats'], transcripts['len_stats'], transcripts['unext_reads'], transcripts['any_reads'] = zip(*transcripts.index.map(read_pos_stats))
(transcripts['mean_5p_ext'],
 transcripts['stdev_5p_ext'],
 transcripts['mode_5p_ext'],
 transcripts['5p_ext_reads'],
 transcripts[f'{count_5p_ext}nt_5p_ext_reads'],
 transcripts['5p_unext_reads'],
 transcripts[f'5p_ext_min_len_{count_5p_ext_min_len}_nt_reads'],
 transcripts[f'{count_5p_ext}nt_5p_ext_min_len_{count_5p_ext_min_len}_nt_reads'],) = zip(*transcripts['5p_ext_stats'])
(transcripts['mean_3p_ext'],
 transcripts['stdev_3p_ext'],
 transcripts['mode_3p_ext'],
 transcripts['3p_ext_reads'],
 transcripts[f'{count_3p_ext}nt_3p_ext_reads'],
 transcripts['3p_unext_reads'],
 transcripts[f'3p_ext_min_len_{count_3p_ext_min_len}_nt_reads'],
 transcripts[f'{count_3p_ext}nt_3p_ext_min_len_{count_3p_ext_min_len}_nt_reads'],) = zip(*transcripts['3p_ext_stats'])
(transcripts['mean_len'],
 transcripts['stdev_len'],
 transcripts['mode_len']) = zip(*transcripts['len_stats'])

del transcripts['read_distrib']
del transcripts['5p_ext_stats']
del transcripts['3p_ext_stats']
del transcripts['len_stats']

transcripts.loc[:, transcripts.dtypes.ne('category')] = transcripts.select_dtypes(exclude='category').fillna(0.0)

transcripts.to_csv(os.path.join(prefix, f"{sample}_transcripts_extension.csv"))
