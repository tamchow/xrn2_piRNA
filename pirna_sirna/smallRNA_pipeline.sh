#!/usr/bin/env bash

export ext='fq.gz'
export min_read_len=16
export max_read_len=40
export min_base_quality=15
export max_unqualified_base_percent=5
export mismatches_allowed=1
export threads=$(($(nproc) / 2))
export genome_fa='c_elegans.PRJNA13758.WS280.genomic.fa'
export index_prefix=$1
export input_prefix=$2
export sample_prefix=$3
export output_prefix=$4

export spike_info="$output_prefix/spike_info.txt"

# external spike is hsa-miR-122
export spike_seq='TGGAGTGTGACAATGGTGTTTG'

export samples=("${@:5}")
export use_parallel=false
export process_siRNAs=true
export rna_mode='mirna'
export zipper="gzip"
if [ ! "$(command -v pigz >/dev/null 2>&1)" ]; then
	# return code != 0 means pigz does not exist
	export zipper="pigz"
fi

if [ -f "$index_prefix.1.ebwt" ] && grep -Fxq $spike_seq $genome_fa; then
	echo "Genome index with spike exists, skipping bowtie-build (delete $index_prefix.*.ebwt to make the index again)."
else
	echo ">Spike" >> $genome_fa
	echo "$spike_seq" >> $genome_fa
	bowtie-build "$genome_fa" "$index_prefix"
fi

pipeline() {
	local i
	i=$1
	local sample
	sample="$sample_prefix$i"
	local infile
	infile="$input_prefix$sample.$ext"
	
	echo "Running adapter detection and read trimming (min length $min_read_len) on $sample"
	
	if [ -n "$ADAPTER_SEQ" ]; then
		fastp --in1 "$infile" --stdout -a "$ADAPTER_SEQ" -q "$min_base_quality" -u "$max_unqualified_base_percent" --length_required "$min_read_len" --length_limit "$max_read_len" --thread "$threads" |
		bowtie -p "$threads" -n "$mismatches_allowed" --tryhard --all --best --strata -S --no-unal "$index_prefix" - |
		samtools view -b -@ $threads |
		samtools sort -l 9 -@ "$threads" > "$output_prefix$sample.sorted.bam"
	elif [ -n "$ADAPTER_FASTA" ]; then
		fastp --in1 "$infile" --stdout --adapter_fasta "$ADAPTER_FASTA" -q "$min_base_quality" -u "$max_unqualified_base_percent" --length_required "$min_read_len" --length_limit "$max_read_len" --thread "$threads" |
		bowtie -p "$threads" -n "$mismatches_allowed" --tryhard --all --best --strata -S --no-unal "$index_prefix" - |
		samtools view -b -@ $threads |
		samtools sort -l 9 -@ "$threads" > "$output_prefix$sample.sorted.bam"
	else
		fastp --in1 "$infile" --stdout -q "$min_base_quality" -u "$max_unqualified_base_percent" --length_required "$min_read_len" --length_limit "$max_read_len" --thread "$threads" |
		bowtie -p "$threads" -n "$mismatches_allowed" --tryhard --all --best --strata -S --no-unal "$index_prefix" - |
		samtools view -b -@ $threads |
		samtools sort -l 9 -@ "$threads" > "$output_prefix$sample.sorted.bam"
	fi
	
	samtools index -@ "$threads" "$output_prefix$sample.sorted.bam"

	local total_reads
	total_reads=$(samtools view -F 260 -@ $threads "$output_prefix$sample.sorted.bam" | cut -f1 | sort -u | wc -l)
	local total_alns
	total_alns=$(samtools view -c -F 260 -@ $threads "$output_prefix$sample.sorted.bam")

	echo "$total_reads mapped unique reads (of $total_alns total alignments) in $sample" | tee -a "$spike_info"
	
	echo "Read alignment, sorting, and indexing completed for $sample, now doing read assignment"

	python ./readext.py "$sample" "$output_prefix" | "$zipper" - > "$output_prefix$sample.assigned_reads.txt.gz"

	if $process_siRNAs; then
		echo "Read assignment completed for $sample, now quantifying siRNA reads"
		python ./quant_siRNAs.py "$sample" "$output_prefix" | "$zipper" - > "$output_prefix$sample.siRNA_reads.txt.gz"
	else
		echo "Read assignment completed for $sample"
	fi;

	echo "Moving on to next sample"
	echo ""
}

len_samples="${#samples[@]}"

if $use_parallel; then
	SHELL="$(type -p bash)"
	export SHELL
	export threads=$((threads / len_samples))
	export threads=$((threads < 1 ? 1 : threads))
	export -f pipeline
	echo "Preparing to process $len_samples samples in parallel using $threads threads per job"
	parallel 'pipeline {}' ::: "${samples[@]}"
else
	echo "Preparing to process $len_samples samples sequentially using $threads threads per job"
	for i in "${samples[@]}"; do
		pipeline "$i"
	done;
fi;

echo "Preprocessing and alignment completed for all samples"