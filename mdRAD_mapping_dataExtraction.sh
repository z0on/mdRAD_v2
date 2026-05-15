
#---- mapping mdRAD data (make sure to filter them first)

REF="/path/to/mygenome.fasta"
>maps
for F in `ls *fastq`; do
echo "bowtie2 --no-unal -x $REF -U $F -S ${F/.fastq/}.sam && samtools sort -O bam -o ${F/.fastq/}.bam ${F/.fastq/}.sam && samtools index ${F/.fastq/}.bam">>maps;
done
# execute all lines in maps. For TACC LAUNCHER system:
ls6_launcher_creator.py -j maps -n maps -t 2:00:00 -w 12 -a IBN21018
sbatch maps.slurm

# extracting methylation counts: collecting reads with identifiable CpG that has been cut my MspJI and writing bismark's cov file for them:
# must have samtools and bedtools

REF="/path/to/mygenome.fasta"

for BAM in `ls *bam`; do
# forward reads:
samtools view $BAM |awk -v OFS='\t' '{if ($2==0 && $5>=20) print $3,$4-14,$4+16}'  | bedtools getfasta -fi $REF -bed - -tab | awk -F'[:-]' '{ print $1,$2,$3}' | awk '{ if (($4~/^CG[GATC][GA]/) && ($4~/[CT][GATC]CG$/)) print $1,$2,0; else if ($4~/^CG[GATC][GA]/) print $1,$2+1,1; else if ($4~/[CT][GATC]CG$/ ) print $1,$3,2 }'  >${BAM}.fwd;
# reverse reads:
samtools view $BAM | awk -v OFS='\t'  '{if ($2==16 && $5>=20) print $3,$4+length($10)-18,$4+length($10)+12}' | bedtools getfasta -fi $REF -bed - -tab | awk -F'[:-]' '{ print $1,$2,$3}' | awk '{ if (($4~/^CG[GATC][GA]/) && ($4~/[CT][GATC]CG$/)) print $1,$2+1,0; else if ($4~/^CG[GATC][GA]/) print $1,$2+1,1; else if ($4~/[CT][GATC]CG$/ ) print $1,$3,2 }'  >${BAM}.rev;
# combining directs and reverses:
cat ${BAM}.fwd ${BAM}.rev | sort | uniq -c | awk -v OFS='\t' '{print $2,$3,$3,100,$1,$4}' | gzip >${BAM}.cov.gz;
rm ${BAM}.fwd;
rm ${BAM}.rev;
done

# this loop works pretty quickly, surprisingly enough, so just run it in idev

# The resulting cov file can be viewed in SeqMonk genome browser. The file has two colums 
# for methylation counts, which are interpreted by SeqMonk as forward and reverse 
# reads one base long. The first column ("forward reads") is actual count of reads 
# for the specific CpG. The second column has values identifying which one of the two 
# MspJI-MseI fragments is being sequenced: 1 (right frag) and 2 (left frag). Value 
# of 0 means that there are two CpGs (both left and right) that might be responsible.
