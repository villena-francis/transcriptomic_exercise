cd part_1

# FastQC (raw data)
mkdir -p quality_controls/fastqc_reports

fastq_dir="data_raw"
fastqc_dir="quality_controls/fastqc_reports"

for fastq in $fastq_dir/*fastq
do
    fastq_base=$(basename "$fastq" .fastq)

    if [ -f "$fastqc_dir/${fastq_base}_fastqc.html" ]
    then
        echo -e "${blue}FastQC analysis already done for $fastq${endblue}"
    else
        fastqc $fastq -o $fastqc_dir
    fi
done

#FastQScreen 
# mkdir -p quality_controls/fastq_screen_reports

# fastq_screen_dir="quality_controls/fastq_screen_reports"

# if [ -f genomes/FastQ_Screen_Genomes/fastq_screen.conf ]
# then
#     echo -e "${blue}Ref genomes already downloaded${endblue}"   
# else
#     fastq_screen --get_genomes --outdir genomes/
# fi

find $fastq_dir -name "*.fastq" | grep -oP 'SRR\d+' | sort | uniq > $fastq_dir/sid_list.txt 

list=$(cat $fastq_dir/sid_list.txt)

# for sid in $list
# do
#     line12=$(sed -n '12p' genomes/FastQ_Screen_Genomes/fastq_screen.conf)
#     if [[ $line12 =~ ^#.* ]]
#     then
#         echo
#         echo -e "${red}Set bowtie2 path in fastq_screen.conf${endred}"
#         echo -e "${red}and then re-run this pipeline${endred}"
#         exit 1
#     else 
#         fastq_screen --conf genomes/FastQ_Screen_Genomes/fastq_screen.conf \
#                      --outdir $fastq_screen_dir \
#                      --aligner bowtie2 \
#                      $fastq_dir/${sid}.chr21_1.fastq $fastq_dir/${sid}.chr21_2.fastq
#     fi
# done

# Cutadapt
mkdir -p cutadapt

trimmed_dir="cutadapt"

for sid in $list
do
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -m 50 -u 15 -U 15 \
         -g CTTTTACTTCCTCTAGATAGTCAAGTTCGACCGTCTTCTCAGCGCTCCGC \
         -o "$trimmed_dir/${sid}.chr21_1_trimmed.fastq" \
         -p "$trimmed_dir/${sid}.chr21_2_trimmed.fastq" \
         $fastq_dir/${sid}.chr21_1.fastq $fastq_dir/${sid}.chr21_2.fastq \
         &> $trimmed_dir/cutadapt.stats
done


# FastQC (trimmed data)
for fastq in $trimmed_dir/*fastq
do
    fastq_base=$(basename "$fastq" .fastq)

    if [ -f "$fastqc_dir/${fastq_base}_fastqc.html" ]
    then
        echo -e "${blue}FastQC analysis already done for $fastq${endblue}"
    else
        fastqc $fastq -o $fastqc_dir 
    fi
done


# GENOME INDEX
if [ -f genomes/hg38_chr21.8.ht2 ]
then
    echo -e "${blue}Genome index already exists${endblue}"
else
    hisat2-build --seed 123 -p 10 genomes/*chromosome.21.fa genomes/hg38_chr21 \
    &> genomes/hisat2_index.log
fi


# READ MAPPING AND QUANTIFICATION
mkdir -p hisat2
mkdir -p htseq
aligments_dir="hisat2"
quant_dir="htseq"

for sid in $list
do
    hisat2 --new-summary --summary-file ${aligments_dir}/hisat2_${sid}.stats \
           --seed 123 --phred33 -p 10 -k 1 \
            -x genomes/hg38_chr21 \
            -1 ${trimmed_dir}/${sid}.chr21_1_trimmed.fastq \
            -2 ${trimmed_dir}/${sid}.chr21_2_trimmed.fastq \
            -S ${aligments_dir}/${sid}_trimmed.sam 
    
    samtools view -bS ${aligments_dir}/${sid}_trimmed.sam > ${aligments_dir}/${sid}_trimmed.bam
    samtools sort ${aligments_dir}/${sid}_trimmed.bam -o ${aligments_dir}/${sid}_trimmed_sorted.bam
    samtools index ${aligments_dir}/${sid}_trimmed_sorted.bam

    htseq-count --format=bam \
                --stranded=reverse \
                --mode=intersection-nonempty \
                --minaqual=10 \
                --type=exon \
                --idattr=gene_id --additional-attr=gene_name \
                ${aligments_dir}/${sid}_trimmed_sorted.bam \
                genomes/Homo_sapiens.GRCh38.109.chr21.gtf > ${quant_dir}/${sid}.htseq
done