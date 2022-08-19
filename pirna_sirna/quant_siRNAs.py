import pandas as pd
import pysam

import os
import sys

sample = sys.argv[1]
prefix = sys.argv[2] if len(sys.argv) > 2 else '.'
threads = int(os.environ.get("threads", len(os.sched_getaffinity(0)) // 2))
siRNA_mode = os.environ.get("sirna_mode", "22g")

sample_smallRNA_reads = pd.read_csv(f'{prefix}/{sample}.assigned_reads.txt.gz', sep=' ', header=None,
                                    names=('name', 'aln_start', 'aln_end', 'aln_len', 'query_len'))

sample_smallRNA_reads = set(sample_smallRNA_reads['name'])

transcripts = pd.read_csv("ref/all_genes.tsv", sep='\t',
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
# sRNA_transcripts = pd.read_csv("mirnas_pirnas.tsv", sep='\t',
#                           index_col='transcript_id',
#                           dtype={'transcript_id': str,
#                                  'Chr': 'category',
#                                  'Start': int,
#                                  'End': int,
#                                  'Strand': 'category',
#                                  'Length': int,
#                                  'gene_id': str,
#                                  'biotype': 'category',
#                                  'name': str,
#                                  })
# transcripts.append(sRNA_transcripts)
transcripts = transcripts[~transcripts.index.duplicated(keep='first')]

transcripts['Start'] = transcripts['Start'] - 1
transcripts['End'] = transcripts['End']

if 'name' not in transcripts.columns:
    transcripts['name'] = transcripts.index

sam_file = pysam.AlignmentFile(os.path.join(prefix, f"{sample}.sorted.bam"), 'rb', threads=threads)


def strand_mismatch(template_strand, query):
    query_strand = '-' if query.is_reverse else '+'
    return query_strand != template_strand


def transcript_22g_siRNAs(transcript_id):
    transcript_data = transcripts.loc[transcript_id]
    reads_22g, reads_26g = 0, 0
    try:
        # print(transcript_id, *transcript_data.to_list(), file=sys.stderr)
        start, end, strand = transcript_data['Start'], transcript_data['End'], transcript_data['Strand']
        for read in sam_file.fetch(transcript_data['Chr'], start, end):
            if read.query_name not in sample_smallRNA_reads:
                if read.get_forward_sequence()[read.qstart: read.qend][0] == 'G' and strand_mismatch(strand, read):
                    aln_length = read.query_alignment_length
                    if 20 <= aln_length <= 23:
                        reads_22g += 1
                    elif 24 <= aln_length <= 27:
                        reads_26g += 1
                    print(transcript_data['name'], read.query_name, reads_22g, reads_26g)
    except Exception as e:
        print(f'Error {str(e)} processing transcript_id {transcript_id}', file=sys.stderr)
        raise e
    return reads_22g, reads_26g


def transcript_all_siRNAs(transcript_id):
    transcript_data = transcripts.loc[transcript_id]
    reads = 0
    try:
        # print(transcript_id, *transcript_data.to_list(), file=sys.stderr)
        start, end, strand = transcript_data['Start'], transcript_data['End'], transcript_data['Strand']
        for read in sam_file.fetch(transcript_data['Chr'], start, end):
            if read.query_name not in sample_smallRNA_reads and strand_mismatch(strand, read):
                aln_length = read.query_alignment_length
                if 20 <= aln_length <= 30:
                    reads += 1
                print(transcript_data['name'], read.query_name, reads)
    except Exception as e:
        print(f'Error {str(e)} processing transcript_id {transcript_id}', file=sys.stderr)
        raise e
    return reads


if siRNA_mode == '22g':
    reads_22g, reads_26g = zip(*transcripts.index.map(transcript_22g_siRNAs))

    transcripts['22G_siRNA_reads'] = reads_22g
    transcripts['26G_siRNA_reads'] = reads_26g
else:
    reads_all = transcripts.index.map(transcript_all_siRNAs)

    transcripts['siRNA_reads'] = reads_all

transcripts.to_csv(os.path.join(prefix, f"{sample}_transcripts_siRNAs_{siRNA_mode}.csv"))
