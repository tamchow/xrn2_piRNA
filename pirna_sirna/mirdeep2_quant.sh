#!/usr/bin/env bash

prefix=$1;
min_read_len=16;
max_read_len=40;
min_base_quality=15;
max_unqualified_base_percent=5;
tmpdir="/tmp/mirdeep2/";
ext="${EXTENSION:-fq}"

ref_hairpin=$2;
ref_mature=$3;
species="${SPECIES:-hsa}";
threads=$(($(nproc) / 2))
samples=("${@:4}");
base_prefix=$(basename "$prefix");
len_samples="${#samples[@]}";
echo "About to process $len_samples samples";

join_by() { 
	local d=${1-} f=${2-};
	if shift 2; then 
		printf %s "$f" "${@/#/$d}";
	fi;
}

all_samples=$(join_by ',' "${samples[@]}");
output_file="run_results_$base_prefix$all_samples.csv";

if [ -f "$tmpdir" ]; then
	rm -r "$tmpdir";
else
	mkdir -p "$tmpdir";
fi;

if [[ "$ref_hairpin" != *"$species"* ]]; then
	echo "Processing hairpin sequences for species $species"
	fname_ref_hairpin=$(basename -s ".fa" "$ref_hairpin")
	extract_mirnas.pl "$ref_hairpin" "$species" > "$tmpdir${fname_ref_hairpin}_$species.fa"
	ref_hairpin="$tmpdir$fname_ref_hairpin"
fi
if [[ "$ref_mature" != *"$species"* ]]; then
	echo "Processing mature sequences for species $species"
	fname_ref_mature=$(basename -s ".fa" "$ref_mature")
	extract_mirnas.pl "$ref_mature" "$species" > "$tmpdir${fname_ref_mature}_$species.fa"
	ref_mature="$tmpdir$fname_ref_mature"
fi

for i in "${samples[@]}"; do
	sample="$prefix$i.$ext.gz";
	tmp_file="$tmpdir$base_prefix$i";

	echo "Starting miRNA quantification for sample $sample";

	if [ -n "$ADAPTER_SEQ" ]; then
		fastp --in1 "$sample" --out1 "$tmp_file.trim.fastq.gz" -a "$ADAPTER_SEQ" -q "$min_base_quality" -u "$max_unqualified_base_percent" --length_required "$min_read_len" --length_limit "$max_read_len" --thread "$threads" --dont_eval_duplication;
	elif [ -n "$ADAPTER_FASTA" ]; then
		fastp --in1 "$sample" --out1 "$tmp_file.trim.fastq.gz" --adapter_fasta "$ADAPTER_FASTA" -q "$min_base_quality" -u "$max_unqualified_base_percent" --length_required "$min_read_len" --length_limit "$max_read_len" --thread "$threads" --dont_eval_duplication;
	else
		fastp --in1 "$sample" --out1 "$tmp_file.trim.fastq.gz" -q "$min_base_quality" -u "$max_unqualified_base_percent" --length_required "$min_read_len" --length_limit "$max_read_len" --thread "$threads" --dont_eval_duplication;
	fi
	if [ "$(command -v pigz >/dev/null 2>&1)" ]; then
		# return code != 0 means pigz does not exist, fallback to normal gzip
		gzip -d -c "$tmp_file.trim.fastq.gz" > "$tmp_file.trim.fastq";
	else
		# pigz exists
		pigz -d -c "$tmp_file.trim.fastq.gz" > "$tmp_file.trim.fastq";
	fi;
	mapper.pl "$tmp_file.trim.fastq" -e -m -l "$min_read_len" -h -s "$tmp_file.collapsed.fa" -v;
	# head -n 6 "$tmp_file.collapsed.fa"
	quantifier.pl -d -j -W -p "$ref_hairpin" -m "$ref_mature" -P -r "$tmp_file.collapsed.fa" -y "$base_prefix$i";
	python ./combine_miRDeep2_results.py "./miRNAs_expressed_all_samples_$base_prefix$i.csv" "$base_prefix$i" "$output_file";
	rm -r "./mapper_logs";
	rm -r "./expression_analyses/expression_analyses_$base_prefix$i/";
	echo "Done quantifying miRNAs in $sample";
done;

echo "Combined quantification of all samples stored in $output_file";