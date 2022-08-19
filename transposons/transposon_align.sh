#!/usr/bin/env bash

hiseq_dist_heuristic=100
# novaseq_dist_heuristic is a config option, may be used later
# shellcheck disable=SC2034
novaseq_dist_heuristic=2000
dist_heuristic=$hiseq_dist_heuristic
genome_index_dir=$1
transposon_index_dir=$2
prefix=$3
annotation=$4
samples=("${@:5}")
threads=$(($(nproc) / 2))
fracOverlap=0.01

# Run from same directory
# ./transposon_align.sh ../../ref/index/cel.genome.ws280 ../../ref/index/cel.transposons.ws280 ../../../../data/cel/lifecycle/gDNA/ ../../ref/cel.geneset_transposons.ws280.gtf.gz ../../../../data/cel/lifecycle/gDNA/

if [[ "${#samples[@]}" == 1 && -d "${samples[0]}" ]]; then
	echo "Got a folder for samples, checking contents..."
	for temp_file in "${prefix}"/*.to_transposons.*; do
		rm "$temp_file"
	done
	for temp_file in "${prefix}"/*.to_genome.*; do
		rm "$temp_file"
	done
	# word splitting is desired here
	# shellcheck disable=SC2207
	samples=($(basename -a "${samples[0]}"/* | grep "_R[12].fq.gz" | sed s/"_R[12].fq.gz"// | sort -u))
	echo "Found samples: ${samples[*]}"
fi

for sample in "${samples[@]}"; do
	if [[ -f "${prefix}/${sample}_R1.fq.gz" && -f "${prefix}/${sample}_R2.fq.gz" ]]; then
		echo "Raw FASTQ.GZ files exist, adapter trimming, base correcting with fastp"
		fastp --in1 "${prefix}/${sample}_R1.fq.gz" --in2 "${prefix}/${sample}_R2.fq.gz" \
			--out1 "${prefix}/${sample}_R1.trim.fq.gz" --out2 "${prefix}/${sample}_R2.trim.fq.gz" \
			-q 20 -u 10 --detect_adapter_for_pe --dont_eval_duplication --correction --thread $threads

		aligner='bowtie2'

		echo "Aligning with $aligner"

		transposon_aligned_file="${prefix}/${sample}.to_transposons.bam"
		transposon_aligned_reads="${prefix}/${sample}.to_transposons._R%.fq.gz"

		$aligner -x "${transposon_index_dir}" -1 "${prefix}/${sample}_R1.trim.fq.gz" -2 "${prefix}/${sample}_R2.trim.fq.gz" --very-sensitive-local --dovetail --al-conc-gz "$transposon_aligned_reads" -p $threads |
			samtools collate -f -@ $threads -O - |
			samtools fixmate -rcm -@ $threads - - |
			samtools sort -l 9 -@ $threads |
			samtools markdup -l 150 -d "${dist_heuristic}" -rS -@ $threads - - >"$transposon_aligned_file"
		samtools index -@ $threads "$transposon_aligned_file"

		genome_aligned_file="${prefix}/${sample}.to_genome.bam"
		genome_aligned_reads="${prefix}/${sample}.to_genome._R%.fq.gz"

		$aligner -x "${genome_index_dir}" -1 "${prefix}/${sample}.to_transposons._R1.fq.gz" -2 "${prefix}/${sample}.to_transposons._R2.fq.gz" --very-sensitive --dovetail --al-conc-gz "$genome_aligned_reads" -p $threads |
			samtools collate -f -@ $threads -O - |
			samtools fixmate -rcm -@ $threads - - |
			samtools sort -l 9 -@ $threads |
			samtools markdup -l 150 -d "${dist_heuristic}" -rS -@ $threads - - >"$genome_aligned_file"
		samtools index -@ $threads "$genome_aligned_file"

		$aligner -x "${transposon_index_dir}" -1 "${prefix}/${sample}.to_genome._R1.fq.gz" -2 "${prefix}/${sample}.to_genome._R2.fq.gz" --very-sensitive-local --dovetail -p $threads |
			samtools collate -f -@ $threads -O - |
			samtools fixmate -rcm -@ $threads - - |
			samtools sort -l 9 -@ $threads |
			samtools markdup -l 150 -d "${dist_heuristic}" -rS -@ $threads - - >"$transposon_aligned_file"
		samtools index -@ $threads "$transposon_aligned_file"

		featureCounts -T $threads -a "$annotation" -p --countReadPairs -C -M -O --fracOverlap $fracOverlap --fraction --ignoreDup --largestOverlap -o "$genome_aligned_file.tsv" -R "CORE" "$genome_aligned_file"

		python assign_transposons.py "$annotation" "$genome_aligned_file.featureCounts" "$genome_aligned_file" "$transposon_aligned_file" $fracOverlap "${prefix}/${sample}.transposon_info.csv"
	else
		echo "Unsupported file found, skipping..."
	fi
done
