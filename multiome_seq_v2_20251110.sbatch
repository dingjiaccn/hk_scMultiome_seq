#!/bin/bash

#SBATCH --job-name=run_multiome_processing_v2
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=cusanovichlab
#SBATCH -o %x.out
#SBATCH --mail-type=END

module load samtools/1.10
module load cutadapt/4.7
module load bowtie2/2.5.3
module load bedtools
module load pigz

export PATH=/home/u15/jiachengd/anaconda3/envs/python39/bin:$PATH #### modify here

home_dir="/xdisk/darrenc/jiachengd/20250625_multiome_seq" #### modify here
output_dir="${home_dir}/run_multiome_processing_v2_20251110" #### modify here
bcl_dir=/xdisk/darrenc/hjkang/MultiomeSeq_111025/BCL #### modify here

# dir not exist (output dir)
fastq_dir="${output_dir}/fastqs"
add_bc_dir="${output_dir}/add_bc"
correct_bc_dir="${output_dir}/correct_bc"
trimmomatic_dir="${output_dir}/trim_fastq"
map_dir="${output_dir}/map_distance"
dedup_dir="${output_dir}/dedup_distance"
counts_dir="${output_dir}/counts_table"
peak_dir="${output_dir}/peak"

mkdir -p "${output_dir}"
mkdir -p "${fastq_dir}"
mkdir -p "${add_bc_dir}"
mkdir -p "${correct_bc_dir}"
mkdir -p "${trimmomatic_dir}"
mkdir -p "${map_dir}"
mkdir -p "${dedup_dir}"
mkdir -p "${counts_dir}"
mkdir -p "${peak_dir}"

add_bc_py="${home_dir}/python_scripts/add_bc_v2_20251110.py"
correct_bc_py="${home_dir}/python_scripts/correct_bc_v2_20251110.py"
trimmomatic="/home/u15/jiachengd/Trimmomatic-0.36/trimmomatic-0.36.jar" #### modify here
trimmomatic_adaptor="/home/u15/jiachengd/Trimmomatic-0.36/adapters" 
java="/home/u15/jiachengd/jdk-21.0.1/bin/java" #### modify here 
star="/groups/darrenc/Jiacheng/cellranger-8.0.1/lib/bin/STAR" 
add_tag_py="${home_dir}/python_scripts/add_tag_v2_20251110.py"
atac_dedup_py="${home_dir}/python_scripts/sc_atac_true_dedup_v2_20251110.py"
featurecounts="/home/u15/jiachengd/subread-2.0.8-Linux-x86_64/bin/featureCounts" #### modify here
macs2="/home/u15/jiachengd/anaconda3/envs/python39/bin/macs2" #### modify here
sinto="/home/u15/jiachengd/anaconda3/envs/python39/bin/sinto"
bcl2fastq="/home/u15/jiachengd/anaconda3/envs/python39/bin/bcl2fastq"

bcl2fastq_samplesheet="${home_dir}/sbatch_jobs/samplesheet.csv"
bc_correct_rna_samplesheet="${home_dir}/python_scripts/bc_correct_rna_samplesheet.txt"
bc_correct_atac_samplesheet="${home_dir}/python_scripts/bc_correct_atac_samplesheet.txt"
wellbc_whitelist="${home_dir}/python_scripts/well_bc_whitelist.txt"
samplebc_whitelist="${home_dir}/python_scripts/sample_bc_whitelist.txt"
cellbcwhitelist="${home_dir}/python_scripts/737K-cratac-v1.txt"

rna_raw_fastq_prefix="cDNA" #### modify here
rna_raw_fastq_flow="S2" #### modify here
atac_raw_fastq_prefix="ATAC" #### modify here
atac_raw_fastq_flow="S1" #### modify here

rna_base="RNA"
atac_base="ATAC"

# ====== run bcl2fastq ========
sed -i 's/<Read Number="3" NumCycles="30" IsIndexedRead="Y"\/>/<Read Number="3" NumCycles="30" IsIndexedRead="N"\/>/g' ${bcl_dir}/RunInfo.xml #### modify here
#### modify basemask as needed
${bcl2fastq} \
    --use-bases-mask=Y65n*,I8n*,Y16n*,Y65n* \
    --create-fastq-for-index-reads \
    --minimum-trimmed-read-length=0 \
    --mask-short-adapter-reads=0 \
    --ignore-missing-positions \
    --ignore-missing-filter \
    --ignore-missing-bcls \
    --no-lane-splitting \
    -R ${bcl_dir} \
    --output-dir=${fastq_dir} \
    --interop-dir=${fastq_dir} \
    --sample-sheet=${bcl2fastq_samplesheet} \
    > ${fastq_dir}/bcl2fastq.log 2>&1
sed -i 's/<Read Number="3" NumCycles="30" IsIndexedRead="N"\/>/<Read Number="3" NumCycles="30" IsIndexedRead="Y"\/>/g' ${bcl_dir}/RunInfo.xml #### modify here
cd ${fastq_dir}
for file in *_R3_*.fastq.gz; do
    mv "$file" "${file/_R3_/_TMP_}"
done
# Then, rename all R2 files to I2
for file in *_R2_*.fastq.gz; do
    mv "$file" "${file/_R2_/_I2_}"
done
# Finally, rename the temporary files (originally R3) to R2
for file in *_TMP_*.fastq.gz; do
    mv "$file" "${file/_TMP_/_R2_}"
done

# =========== append bc ==========
echo 'formating RNA fastq...'
python "${add_bc_py}" \
  -1 <(zcat "${fastq_dir}/${rna_raw_fastq_prefix}_${rna_raw_fastq_flow}_R1_001.fastq.gz") \
  -2 <(zcat "${fastq_dir}/${rna_raw_fastq_prefix}_${rna_raw_fastq_flow}_R2_001.fastq.gz") \
  -I1 <(zcat "${fastq_dir}/${rna_raw_fastq_prefix}_${rna_raw_fastq_flow}_I1_001.fastq.gz") \
  -I2 <(zcat "${fastq_dir}/${rna_raw_fastq_prefix}_${rna_raw_fastq_flow}_I2_001.fastq.gz") \
  -a RNA \
  -o "${add_bc_dir}/${rna_base}"
pigz -p 10 "${add_bc_dir}/${rna_base}_R1.fastq"
pigz -p 10 "${add_bc_dir}/${rna_base}_R2.fastq"
echo 'done formating RNA fastq'

echo 'formating ATAC fastq...'
python "${add_bc_py}" \
  -1 <(zcat "${fastq_dir}/${atac_raw_fastq_prefix}_${atac_raw_fastq_flow}_R1_001.fastq.gz") \
  -2 <(zcat "${fastq_dir}/${atac_raw_fastq_prefix}_${atac_raw_fastq_flow}_R2_001.fastq.gz") \
  -I1 <(zcat "${fastq_dir}/${atac_raw_fastq_prefix}_${atac_raw_fastq_flow}_I1_001.fastq.gz") \
  -I2 <(zcat "${fastq_dir}/${atac_raw_fastq_prefix}_${atac_raw_fastq_flow}_I2_001.fastq.gz") \
  -a ATAC \
  -o "${add_bc_dir}/${atac_base}"
pigz -p 10 "${add_bc_dir}/${atac_base}_R1.fastq"
pigz -p 10 "${add_bc_dir}/${atac_base}_R2.fastq"
echo 'done formating ATAC fastq'

# =========== correct bc ==========
echo 'correcting RNA barcode...'
python "${correct_bc_py}" \
  -1 <(zcat "${add_bc_dir}/${rna_base}_R1.fastq.gz") \
  -2 <(zcat "${add_bc_dir}/${rna_base}_R2.fastq.gz") \
  -o "${correct_bc_dir}/${rna_base}" \
  -X \
  --samplesheet "${bc_correct_rna_samplesheet}" \
  --stats_out "${correct_bc_dir}/${rna_base}_bc.fixed.stats" \
  --wellbc "${wellbc_whitelist}" \
  --samplebc "${samplebc_whitelist}" \
  --cellbc "${cellbcwhitelist}" \
  --assay RNA
mv "${correct_bc_dir}/${rna_base}.${rna_base}_R1.fastq" "${correct_bc_dir}/${rna_base}_R1.fastq"
mv "${correct_bc_dir}/${rna_base}.${rna_base}_R2.fastq" "${correct_bc_dir}/${rna_base}_R2.fastq"
pigz -p 10 "${correct_bc_dir}/${rna_base}_R1.fastq"
pigz -p 10 "${correct_bc_dir}/${rna_base}_R2.fastq"

echo 'correcting ATAC barcode...'
python "${correct_bc_py}" \
  -1 <(zcat "${add_bc_dir}/${atac_base}_R1.fastq.gz") \
  -2 <(zcat "${add_bc_dir}/${atac_base}_R2.fastq.gz") \
  -o "${correct_bc_dir}/${atac_base}" \
  -X \
  --samplesheet "${bc_correct_atac_samplesheet}" \
  --stats_out "${correct_bc_dir}/${atac_base}_bc.fixed.stats" \
  --wellbc "${wellbc_whitelist}" \
  --samplebc "${samplebc_whitelist}" \
  --cellbc "${cellbcwhitelist}" \
  --assay ATAC
mv "${correct_bc_dir}/${atac_base}.${atac_base}_R1.fastq" "${correct_bc_dir}/${atac_base}_R1.fastq"
mv "${correct_bc_dir}/${atac_base}.${atac_base}_R2.fastq" "${correct_bc_dir}/${atac_base}_R2.fastq"
pigz -p 10 "${correct_bc_dir}/${atac_base}_R1.fastq"
pigz -p 10 "${correct_bc_dir}/${atac_base}_R2.fastq"

# =========== trimming adpators / poly-As =============
echo 'remove polyA for cDNA sequence'
cutadapt -a "A{20}" \
-o "${trimmomatic_dir}/${rna_base}_R2_trim_polyA.fastq.gz" \
"${trimmomatic_dir}/${rna_base}_R2.fastq.gz" \
--minimum-length 30 \
--overlap 10 \
--cores 0 \
> "${trimmomatic_dir}/${rna_base}_trim_polyA_stat.txt" 2>&1

echo 'trimming ATAC (nextera R1 and R2) adaptors...'
${java} -Xmx1G -jar ${trimmomatic} \
PE \
-threads 10 \
${correct_bc_dir}/${atac_base}_R1.fastq.gz ${correct_bc_dir}/${atac_base}_R2.fastq.gz \
${trimmomatic_dir}/${atac_base}_R1_trimmomatic.fastq.gz ${trimmomatic_dir}/${atac_base}_R1_trimmomatic_unpaired.fastq.gz \
${trimmomatic_dir}/${atac_base}_R2_trimmomatic.fastq.gz ${trimmomatic_dir}/${atac_base}_R2_trimmomatic_unpaired.fastq.gz \
ILLUMINACLIP:${adapter_path}/multiome_atac-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 \
2> ${trimmomatic_dir}/${atac_base}_trimmomatic.log;
echo 'done trimming ATAC adaptors'

# ======= map to genome ========
echo 'star mapping rna reads to genome...'
genome_star_mix=/groups/darrenc/references/10x/refdata-gex-GRCh38-and-mm10-2020-A/star
genome_star_mix_fa=/groups/darrenc/references/10x/refdata-gex-GRCh38-and-mm10-2020-A/fasta/genome.fa
mkdir -p "${map_dir}/${rna_base}"
${star} \
     --runThreadN 10 \
     --genomeDir ${genome_star_mix} \
     --readFilesIn "${trimmomatic_dir}/${rna_base}_R2_trim_polyA.fastq.gz" \
     --readFilesCommand zcat \
     --outFileNamePrefix ${map_dir}/${rna_base}/ \
     --outSAMtype BAM SortedByCoordinate \
     --quantTranscriptomeBan Singleend \
     --quantMode GeneCounts \
     --outSAMattributes NH HI AS NM MD \
     --outFilterMultimapNmax 1
mv "${map_dir}/${rna_base}/Aligned.sortedByCoord.out.bam" "${map_dir}/${rna_base}/${rna_base}.bam"
samtools index "${map_dir}/${rna_base}/${rna_base}.bam" 
echo 'done mapping rna reads'

echo 'bowtie2 mapping atac reads to genome...'
genome_bowtie_mix=/groups/darrenc/references/bowtie/hg38_mm10/hg38_mm10
genome_bowtie_hg=/groups/darrenc/references/bowtie/hg38/hg38
genome_bowtie_mm=/groups/darrenc/references/bowtie/mm10/mm10
mkdir -p "${map_dir}/${atac_base}"
bowtie2 -p 8 -X 2000 -3 1 -x ${genome_bowtie_mix} \
-1 ${trimmomatic_dir}/${atac_base}_R1_trimmomatic.fastq.gz \
-2 ${trimmomatic_dir}/${atac_base}_R2_trimmomatic.fastq.gz \
2> ${map_dir}/${atac_base}/${atac_base}.bowtie2.log | samtools view -bS - > "${map_dir}/${atac_base}/${atac_base}.bam"

echo 'filtering out low quality read pairs and unknown chr contigs...'
samtools view -h -f3 -F12 -q10 "${map_dir}/${atac_base}/${atac_base}.bam" \
| grep -v '[0-9]'$'\t'hg38chrM | grep -v '[0-9]'$'\t'mm10chrM \
| grep -v '[0-9]'$'\t'GL | grep -v '[0-9]'$'\t'KI | grep -v '[0-9]'$'\t'JH \
| samtools view -Su - | samtools sort -@ 16 -T ${map_dir}/${atac_base}/${atac_base}.sorttemp \
-o "${map_dir}/${atac_base}/${atac_base}.q10.bam"
samtools index "${map_dir}/${atac_base}/${atac_base}.q10.bam"

# ====== format bam file =======
echo 'formatting RNA.bam...'
python ${add_tag_py} -i "${map_dir}/${rna_base}/${rna_base}.bam" \
-o "${map_dir}/${rna_base}/${rna_base}_tagged.bam" \
-a RNA
samtools index "${map_dir}/${rna_base}/${rna_base}_tagged.bam"
echo 'done formating'

echo 'formatting ATAC.bam...'
python ${add_tag_py} -i "${map_dir}/${atac_base}/${atac_base}.q10.bam" \
-o "${map_dir}/${atac_base}/${atac_base}_tagged.bam" \
-a ATAC
samtools index "${map_dir}/${atac_base}/${atac_base}_tagged.bam"

# ====== deduplication ======
# try deduplication -unique- method
mkdir -p "${dedup_dir}/${rna_base}"
umi_tools dedup \
  -I "${map_dir}/${rna_base}/${rna_base}_tagged.bam" \
  -S "${dedup_dir}/${rna_base}/${rna_base}_dedup.bam" \
  --extract-umi-method=tag \
  --umi-tag=UB \
  --cell-tag=CB \
  --per-cell \
  --method=unique
echo 'done RNA PCR duplicates'

echo 'deduplicating ATAC reads...'
mkdir -p "${dedup_dir}/${atac_base}"
python ${atac_dedup_py} \
"${map_dir}/${atac_base}/${atac_base}_tagged.bam" \
"${dedup_dir}/${atac_base}/${atac_base}_dedup.bam"
samtools sort -o "${dedup_dir}/${atac_base}/${atac_base}_dedup.sorted.bam" \
  "${dedup_dir}/${atac_base}/${atac_base}_dedup.bam"
samtools index "${dedup_dir}/${atac_base}/${atac_base}_dedup.sorted.bam"
echo 'done ATAC reads deduplication'

# ========= RNA count matrix ========
mkdir -p "${counts_dir}/${rna_base}"
hg_gtf="/groups/darrenc/Jiacheng/10x_hg38_genes.gtf"
mx_gtf="/groups/darrenc/references/10x/refdata-gex-GRCh38-and-mm10-2020-A/genes/genes.gtf"
# ---- annotation ------
${featurecounts} -T 10 \
-a ${mx_gtf} \
-t gene -g gene_name \
-R BAM \
-o ${counts_dir}/${rna_base}/${rna_base}_featureCounts \
"${dedup_dir}/${rna_base}/${rna_base}_dedup.bam"
samtools sort -o "${counts_dir}/${rna_base}/${rna_base}_dedup.bam.featureCounts.sorted.bam" \
  "${counts_dir}/${rna_base}/${rna_base}_dedup.bam.featureCounts.bam"
samtools index "${counts_dir}/${rna_base}/${rna_base}_dedup.bam.featureCounts.sorted.bam"
# ----- counts table -----
# umi + position
rna_counter='/xdisk/darrenc/jiachengd/20250625_multiome_seq/python_scripts/generate_counts_table.py'
python ${rna_counter} \
"${counts_dir}/${rna_base}/${rna_base}_dedup.bam.featureCounts.sorted.bam" \
"${counts_dir}/${rna_base}"

# =========== ATAC count matrix =========
mkdir -p "${counts_dir}/${atac_base}"
${sinto} fragments \
  --bam "${dedup_dir}/${atac_base}/${atac_base}_dedup.sorted.bam" \
  --barcodetag 'CB' \
  --fragments "${counts_dir}/${atac_base}/atac_fragments.tsv" \
  --collapse_within \
  --use_chrom .
  -m 0 \ # do not perform any qc
  --min_distance 0 \ # do not perform any qc
  --max_distance 100000 \ # do not perform any qc
cat "${counts_dir}/${atac_base}/atac_fragments.tsv" | \
sort -k1,1 -k2,2n | \
bgzip -c > "${counts_dir}/${atac_base}/atac_fragments.tsv.gz"
tabix -p bed atac_fragments.tsv.gz

# ========= call ATAC peaks =============
mkdir -p ${peak_dir}
bedtools bamtobed \
-i "${dedup_dir}/${atac_base}/${atac_base}_dedup.sorted.bam" > \
"${dedup_dir}/${atac_base}/${atac_base}_dedup.sorted.bed";

${macs2} callpeak \
-t "${dedup_dir}/${atac_base}/${atac_base}_dedup.sorted.bed" \
-f BED -g 4.6e9 --nomodel --keep-dup all --extsize 200 --shift -100 \
--call-summits -n ${peak_dir}/${atac_base}.nodups;

sort -k 8gr,8gr ${peak_dir}/${atac_base}.nodups_peaks.narrowPeak | \
awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | sort -k1,1 -k2,2n -k3,3n | \
gzip -c > ${peak_dir}/${atac_base}.nodups.narrowPeak.gz;

rm ${peak_dir}/${atac_base}.nodups_peaks.narrowPeak;

zcat ${peak_dir}/${atac_base}.nodups.narrowPeak.gz | \
bedtools merge -i stdin > ${peak_dir}/${atac_base}.nodups.narrowPeak_merged.bed






