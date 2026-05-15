#!/bin/bash

# Help message #####
script_name=`basename $0`
help_msg="Extract mdRAD methylation data

   Extracts CpG locations from mdRAD alignments.
   Outputs a .cov.gz file with two colums  for methylation counts.
   The first column is the actual count of reads for the specific CpG. 
   The second column has values identifying which one of the two 
   MspJI-MseI fragments is being sequenced: 1 (right frag) and 2 (left frag). 
   Value of 0 means that there are two CpGs (both left and right) that might 
   be responsible.
   
Usage: ${script_name} bam genome out_dir
       ${script_name} --help
       
       bam: sorted and indexed alignment in bam format
       genome: fasta file the bam was aligned to
       out_dir: output directory (default: same as bam)
       
       --help: print this message
       
"
# Function Definitions #########
get_md_positions() {
  local BAM=$1
  local direc=${2:-fwd} # other option should be rev
  
  # Both of these get the genome position of the read & return a bed file
  #  w/ the chromosome name, the start, and the end positions
  
  if [ $direc = "fwd" ]; then
    local awk_string='{if ($2==0 && $5>=20) print $3,$4-14,$4+16}'  
  elif [ $direc = "rev" ]; then
    local awk_string='{if ($2==16 && $5>=20) print $3,$4+length($10)-18,$4+length($10)+12}'
  else 
    >&2 echo "Error, direction must be 'fwd' or 'rev'; user supplied $direc"
    exit 1
  fi
  # View bam, remove unaligned, and pull out the positions of the MD sites in reads
  # The second awk sequence removes negative start coordinates
  # Fortunately, bedtools will remove sequences beyond the lenght for me
  samtools view -F 0x04 $BAM | awk -v OFS='\t' "$awk_string" | \
  awk -v OFS='\t' '{if ($2>=0) print $0}'
}
pos_to_seq() {
  local REF=${1:-genome/Amil.fasta}
  # Converts md_positions to the rererence sequence; the awk line re-formats to a tsv
  bedtools getfasta -fi $REF -bed - -tab | \
  awk -F'[:-]' '{ print $1,"\t",$2,"\t",$3}' 
}
get_methyl_score() {
  # Filter sequences to only include those starting / ending w/ an MD_rad sequence
  # 0 = ambiguous
  # 1 = left
  # 2 = right?
  # Left & right patterns to match
  local l_pat='($4~/^CG[GATC][GA]/)'
  local r_pat='($4~/[CT][GATC]CG$/)'
  # NOTE: IN the initial version of this, the second column was just $2 for forward; 
  # I think this was a mistake
  # Result templates
  local l_print='print $1,$2+1' 
  local r_print='print $1,$3'
  # Combined queries
  local zero="($l_pat && $r_pat) ${l_print},0"
  local one="$l_pat ${l_print},1"
  local two="$r_pat ${r_print},2"
  # Final awk text
  local awk_txt="{if  $zero; else if $one; else if $two }"
  awk "$awk_txt"
}
combine_fwd_rev() {
  # Combines separarte forward and reverse reads
  # Not 100% sure on the columns
  local fwd=${1}
  local rev=${2}
  cat $fwd $rev | sort | uniq -c | \
    awk -v OFS='\t' '{print $2,$3,$3,100,$1,$4}' 
}
# Run the whole sequence
extract_mdrad() {
  local bam=${1}
  local genome=${2}
  local covgz_dir=${3:-`dirname $bam`}
  local base_bam=`basename $bam`
  local covgz=${covgz_dir}/${base_bam}.cov.gz
  local tmp=/tmp/$base_bam
  
  mkdir -p $covgz_dir
  # Run forward & reverse extractions
  get_md_positions $bam fwd | pos_to_seq $genome | get_methyl_score > ${tmp}_fwd
  get_md_positions $bam rev | pos_to_seq $genome | get_methyl_score > ${tmp}_rev
  
  # Combine & save
  combine_fwd_rev ${tmp}_fwd ${tmp}_rev | gzip > ${covgz}
}

# Validate argument and run some tests ####

# Test to see if bam and fasta files exist
stop_if_no_file() {
  if [ ! -f $1 ]; then
  >&2 echo "Error: file $1 does not exist" 
  exit 1  
  fi
}
stop_if_missing_arg() {
  arg_name=$1
  arg_val=$2
  if [ -z ${arg_val+x} ]; then # true if arg_val is unset
    >&2 echo "Error: Argument ${arg_name} not provided;"
    >&2 echo "Use '${script_name} --help' for usage information."
    exit 1
  fi
}

stop_if_missing_arg "bam" $1
stop_if_missing_arg "genome" $2


# Print the help message
if [[ $1 = "--help" || $1 = "-h"  ]]; then
  echo "$help_msg"
  exit 0 # or should it be return?
fi

bam_file=$1
fasta_file=$2
out_dir=$3

# Validate files exist
stop_if_no_file $bam_file
stop_if_no_file $fasta_file


# Test to see that bedtools is installed
if [ -z `which bedtools` ]; then
  >&2 echo "Error: cannot find bedtools" 
  exit 1
fi
# Run the script #####
make_mdrad $1 $2 $3
