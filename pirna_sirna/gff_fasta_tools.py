#! /usr/bin/env python3

from enum import Enum
import gzip
from pathlib import Path
import sys
from typing import Iterator, NamedTuple
from Bio.SeqIO import parse as seq_parse, write as seq_write, SeqRecord
from Bio.Seq import Seq
import pandas as pd

GenomeDict = dict[str, Seq]


class Strand(Enum):
    FORWARD = "+"
    REVERSE = "-"

    def __str__(self):
        return str(self.value)

    @staticmethod
    def from_str(strand_literal: str):
        for member in Strand:
            if member.value == strand_literal:
                return member
        valid_strand_types = ", ".join(strand_type.value for strand_type in Strand)
        raise ValueError(
            f"{strand_literal} is not a valid strand type ({valid_strand_types})"
        )

    @staticmethod
    def from_bool(is_reverse: bool = False):
        if is_reverse:
            return Strand.REVERSE
        else:
            return Strand.FORWARD


def load_genome(genome_fasta: Path) -> GenomeDict:
    handle = genome_fasta.open("rt")
    if genome_fasta.suffix.lower() == ".gz":
        handle = gzip.open(genome_fasta, "rt")
    return {record.id: record.seq for record in seq_parse(handle, "fasta")}


def get_sequence(
    genome: GenomeDict, chromosome: str, start: int, end: int, strand: Strand
) -> Seq:
    chromosome_sequence = genome[chromosome]
    start = max(start, 1)
    seq: Seq = chromosome_sequence[start - 1 : end]
    if strand == Strand.REVERSE:
        return seq.reverse_complement()
    else:
        return seq


def load_gff_or_gtf(
    gff_or_gtf_file: Path, only_genes_from_gtf: bool = True
) -> pd.DataFrame:
    addtl_info_cols: set[str] = set()
    is_gtf = any("gtf" in suffix.lower() for suffix in gff_or_gtf_file.suffixes)

    def addtl_info_str_to_dict(addtl_info_str: str) -> dict[str, str]:
        return {
            sub_info_split[0]: (
                sub_info_split[1][1:-1] if is_gtf else sub_info_split[1]
            )  # remove quotes around extra info values in last column of GTF files
            for sub_info in (
                addtl_info_str.split("; ") if is_gtf else addtl_info_str.split(";")
            )
            # splitting with '=' for sub-records should always return a list of 2, which is truthy
            if (
                sub_info_split := (
                    sub_info.split(" ", maxsplit=1) if is_gtf else sub_info.split("=")
                ),
                addtl_info_cols.add(sub_info_split[0]),
            )
        }

    gff_or_gtf = pd.read_table(
        gff_or_gtf_file,
        header=None,
        usecols=(0, 1, 2, 3, 4, 6, 8),
        names=("chromosome", "source", "type", "start", "end", "strand", "addtl_info"),
        comment="#",
    )
    if "length" not in gff_or_gtf.columns:
        gff_or_gtf["length"] = gff_or_gtf["end"] - gff_or_gtf["start"] + 1
    gff_or_gtf["strand"] = gff_or_gtf["strand"].apply(Strand.from_str)
    if only_genes_from_gtf:
        gff_or_gtf = gff_or_gtf[gff_or_gtf["type"].str.lower() == "gene"]
    gff_or_gtf["addtl_info"] = gff_or_gtf["addtl_info"].apply(addtl_info_str_to_dict)
    for col in addtl_info_cols:
        gff_or_gtf[col] = gff_or_gtf["addtl_info"].apply(
            lambda info: info.get(col, None)
        )
    
    id_col = "gene_id" if is_gtf else "ID"
    gff_or_gtf[id_col] = normalize_transposon_id(gff_or_gtf, id_col)
    del gff_or_gtf["addtl_info"]
    return gff_or_gtf


def normalize_transposon_id(df: pd.DataFrame, id_col: str) -> pd.Series:
    return df[id_col].apply(
        lambda gene_id: (
            gene_id.replace("Transposon:", "").replace("WBTransposon", "WBT")
        )
    )


def get_sequence_for_gff_record(genome: GenomeDict, gff_record: NamedTuple) -> Seq:
    return get_sequence(
        genome,
        gff_record.chromosome,
        gff_record.start,
        gff_record.end,
        gff_record.strand,
    )


def gff_to_seq_iter(genome: GenomeDict, gff: pd.DataFrame) -> Iterator[SeqRecord]:
    for record in gff.itertuples(index=False):
        description = ",".join(f"{k}={v}" for k, v in record._asdict().items())
        yield SeqRecord(
            get_sequence_for_gff_record(genome, record),
            id=record.ID,
            name=record.Name,
            description=description,
        )


def gff_to_fasta(genome: GenomeDict, gff: pd.DataFrame, fasta_output_file: Path) -> int:
    output_handle = fasta_output_file.open("wt")
    if fasta_output_file.suffix == ".gz":
        output_handle = gzip.open(fasta_output_file, "wt")
    return seq_write(gff_to_seq_iter(genome, gff), output_handle, "fasta")


if __name__ == "__main__":
    gff_to_fasta(
        load_genome(Path(sys.argv[1])),
        load_gff_or_gtf(Path(sys.argv[2])),
        Path(sys.argv[3]),
    )
