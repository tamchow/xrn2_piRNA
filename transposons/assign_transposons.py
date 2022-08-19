from pathlib import Path
import sys, os
from pathlib import Path
from typing import Literal, NamedTuple, Optional

import pandas as pd
import pysam

from gff_fasta_tools import Strand, load_gff_or_gtf, normalize_transposon_id

threads = os.cpu_count()


class GenomicFeature(NamedTuple):
    # transcript_id: Optional[str]
    gene_id: Optional[str]
    gene_biotype: Optional[str]
    chromosome: Optional[str]
    start: Optional[int]
    end: Optional[int]
    length: Optional[int]
    strand: Optional[Strand]

    @staticmethod
    def from_pdSeries(row: pd.Series):
        try:
            return GenomicFeature(*[row[field] for field in GenomicFeature._fields])
        except KeyError as key_error:
            print(row)
            raise key_error

    @staticmethod
    def empty():
        return GenomicFeature(None, None, None, None, None, None, None)


def detail_transposon_feature_assignments(
    gene_info_gtf: Path, featureCounts_read_assignments_tsv: Path
) -> pd.DataFrame:
    gene_info = load_gff_or_gtf(gene_info_gtf, only_genes_from_gtf=True)
    gene_info.set_index("gene_id", inplace=True)
    gene_info["gene_id"] = gene_info.index
    del gene_info["source"]
    del gene_info["type"]
    featureCounts_info = pd.read_table(
        featureCounts_read_assignments_tsv,
        usecols=("read_id", "assigned", "targets"),
        names=("read_id", "assigned", "target_count", "targets"),
    )
    featureCounts_info = featureCounts_info[
        featureCounts_info["assigned"].isin(("Assigned", "Unassigned_NoFeatures"))
    ]

    featureCounts_info["assigned"] = featureCounts_info["assigned"].apply(
        lambda status: status.startswith("Assigned")
    )
    featureCounts_info["targets"] = featureCounts_info[["assigned", "targets"]].apply(
        lambda row: row["targets"] if row["assigned"] else None, axis="columns"
    )

    featureCounts_info["targets"] = featureCounts_info["targets"].apply(
        lambda targets: [
            GenomicFeature.from_pdSeries(gene_info.loc[target])
            for target in targets.split(",")
        ]
        if not pd.isna(targets)
        else []
    )
    max_targets = max(map(len, featureCounts_info["targets"]))
    target_expanded_df = pd.DataFrame(
        featureCounts_info["targets"]
        .apply(
            lambda targets: [
                field
                for target in targets
                + [GenomicFeature.empty()] * (max_targets - len(targets))
                for field in target
            ]  # pad out
        )
        .tolist(),
        columns=[
            f"{field}.{i}"
            for i in range(1, max_targets + 1)
            for field in GenomicFeature._fields
        ],
        index=featureCounts_info.index,
    )
    featureCounts_info = featureCounts_info.join(target_expanded_df, how="inner")
    del featureCounts_info["assigned"]
    del featureCounts_info["targets"]
    return gene_info, featureCounts_info


def detail_transposon_locations(
    reads_to_genome_bam: Path,
    reads_to_transposons_bam: Path,
    gene_info: pd.DataFrame,
    feature_assignments: pd.DataFrame,
) -> pd.DataFrame:
    to_genome = pysam.AlignmentFile(
        str(reads_to_genome_bam), mode="rb", threads=threads
    )
    to_transposons = pysam.AlignmentFile(
        str(reads_to_transposons_bam), mode="rb", threads=threads
    )

    feature_assignments = feature_assignments.set_index("read_id")

    def aln_for_transposon_read(aln):
        return (
            aln.query_name,
            aln.reference_name,
            aln.reference_start,
            aln.reference_end,
            aln.query_alignment_length,
            Strand.from_bool(aln.is_reverse),
            aln.query_alignment_start,
            aln.query_alignment_end,
        )

    to_genome_df = pd.DataFrame(
        data=[aln_for_transposon_read(aln) for aln in to_genome.fetch()],
        columns=[
            "read_id",
            "chromosome",
            "genome_start",
            "genome_end",
            "genome_aln_len",
            "genome_strand",
            "genome_query_start",
            "genome_query_end",
        ],
    )
    to_genome_df = (
        to_genome_df.sort_values(by=["read_id", "genome_aln_len"])
        .drop_duplicates(subset=["read_id"], keep="last")
        .set_index("read_id")
    )

    to_transposons_df = pd.DataFrame(
        data=[
            aln_for_transposon_read(aln) for aln in to_transposons.fetch(until_eof=True)
        ],
        columns=[
            "read_id",
            "transposon",
            "transposon_ref_start",
            "transposon_ref_end",
            "transposon_aln_len",
            "transposon_ref_strand",
            "transposon_query_start",
            "transposon_query_end",
        ],
    )
    to_transposons_df = (
        to_transposons_df.sort_values(by=["read_id", "transposon_aln_len"])
        .drop_duplicates(subset=["read_id"], keep="last")
        .set_index("read_id")
    )
    to_transposons_df["transposon"] = normalize_transposon_id(
        to_transposons_df, "transposon"
    )

    to_transposons_df[
        [
            "transposon_family",
            "transposon_chromosome",
            "transposon_start",
            "transposon_end",
            "transposon_length",
            "transposon_strand",
        ]
    ] = to_transposons_df["transposon"].apply(
        lambda gene_id: gene_info.loc[
            gene_id,
            [
                "gene_name",
                "chromosome",
                "start",
                "end",
                "length",
                "strand",
            ],
        ]
    )

    feature_assignments = feature_assignments.join(to_genome_df, how="inner")
    feature_assignments = feature_assignments.join(to_transposons_df, how="inner")

    def insertion_type(df_row):
        genomic_start, genomic_end, ref_start, ref_end, query_len = (
            df_row["transposon_start"],
            df_row["transposon_end"],
            df_row["genome_start"],
            df_row["genome_end"],
            df_row["genome_aln_len"],
        )

        if not (genomic_start <= ref_start <= ref_end <= genomic_end):
            return "Invalid Coordinates"
        flank_5p = genomic_start, min(genomic_start + query_len, genomic_end)
        flank_3p = max(genomic_start, genomic_end - query_len), genomic_end
        if flank_5p[0] <= ref_start <= flank_5p[1]:
            return "5p"
        if flank_3p[0] <= ref_end <= flank_3p[1]:
            return "3p"
        return "Not Insert"

    feature_assignments["insertion_type"] = feature_assignments.apply(
        insertion_type, axis="columns"
    )

    feature_assignments = feature_assignments[
        feature_assignments["insertion_type"] != "Invalid Coordinates"
    ]
    return feature_assignments


def overlap(ref_start: int, ref_end: int, q_start: int, q_end: int) -> float:
    q_len = q_end - q_start
    if q_start > ref_end or q_end < ref_start:
        return 0
    if q_start <= ref_start <= q_end:
        return -(q_end - ref_start) / q_len
    if q_start <= ref_end <= q_end:
        return (ref_end - q_start) / q_len
    return 1


def collapse_transposon_info(
    transposon_df: pd.DataFrame,
    collapse_col: Literal["transposon", "transposon_family"] = "transposon",
    overlap_frac_threshold_for_collapse: float = 0.5,
    dedup_on_feature_assignments: bool = False,
) -> pd.DataFrame:
    transposon_df_grouped = transposon_df.groupby(collapse_col)
    dedup_df = pd.DataFrame(columns=transposon_df.columns)
    processed = 0
    for key, df in transposon_df_grouped:
        print(f"Have {len(df)} entries for transposon {key}")
        dedup_accum = [df.iloc[0].copy()]
        dedup_accum[0]["collapsed_reads"] = 1
        for _, row in df.iterrows():
            row = row.copy()
            for idx, ref_row in enumerate(dedup_accum):
                if row["chromosome"] == ref_row["chromosome"]:
                    overlap_frac = overlap(
                        ref_row["genome_start"],
                        ref_row["genome_end"],
                        row["genome_start"],
                        row["genome_end"],
                    )
                    if overlap_frac_threshold_for_collapse <= overlap_frac < 1:
                        dedup_accum[idx]["collapsed_reads"] += 1
                        dedup_accum[idx]["genome_end"] = row["genome_end"]
                        break
                    elif -1 < overlap_frac <= -overlap_frac_threshold_for_collapse:
                        dedup_accum[idx]["collapsed_reads"] += 1
                        dedup_accum[idx]["genome_start"] = row["genome_start"]
                        break
                    elif overlap_frac == 1:
                        dedup_accum[idx]["collapsed_reads"] += 1
                        break
            else:
                row["collapsed_reads"] = 1
                dedup_accum.append(row)

        dedup_df = dedup_df.append(dedup_accum)
        processed += 1
        print(
            f"Added {len(dedup_accum)} rows to dedup_df for transposon {key} with {len(df)} initial entries, {processed} entries processed among {len(transposon_df_grouped)}"
        )
    dedup_df.reset_index(drop=True, inplace=True)
    dedup_df.dropna(how="all", axis="columns", inplace=True)
    if dedup_on_feature_assignments:
        dedup_df.drop_duplicates(
            subset=[col for col in dedup_df.columns if "." in col], inplace=True
        )
    dedup_df.drop_duplicates(inplace=True)
    dedup_df["location_category"] = dedup_df[
        [
            "chromosome",
            "genome_start",
            "genome_end",
            "transposon_chromosome",
            "transposon_start",
            "transposon_end",
        ]
    ].apply(
        lambda row: "Expected Location"
        if (row["chromosome"] == row["transposon_chromosome"])
        and (
            abs(
                overlap(
                    row["transposon_start"],
                    row["transposon_end"],
                    row["genome_start"],
                    row["genome_end"],
                )
            )
            >= 0.5
        )
        else "Unexpected Location",
        axis="columns",
    )
    return dedup_df


def test(gene_info, featureCounts, to_genome_bam, to_transposons_bam, overlap_frac_threshold):
    overlap_frac_threshold = float(overlap_frac_threshold)
    gene_info_df, featureCounts_df = detail_transposon_feature_assignments(
        Path(gene_info), Path(featureCounts)
    )
    print("Transposon feature assignments processed")
    feature_assignments_df = detail_transposon_locations(
        Path(to_genome_bam),
        Path(to_transposons_bam),
        gene_info_df,
        featureCounts_df,
    )
    print("Transposon insert locations detailed")
    dedup_df = collapse_transposon_info(feature_assignments_df, overlap_frac_threshold_for_collapse=overlap_frac_threshold)
    print("Deduplication done")
    return dedup_df


def to_clean_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, na_rep="", float_format="%.0f", index=False)


if __name__ == "__main__":
    to_clean_csv(test(*sys.argv[1:6]), sys.argv[6])
