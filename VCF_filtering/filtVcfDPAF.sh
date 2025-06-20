# filtVCF is a function meant to filter VCF files with a
# sequencing depth >= 10, and a population allelic frequency < 0.0001 by default
# different values can be specified in the second and third argument of the function
# the output file is saved as filt_OUTPUT

filtVcfDPAF() {
  module load bcftools
  if [ $#-lt 1 ] || [ $#-gt 3 ]; then
    echo "Usage: filtVcf <input_vcf_file> [min_depth] [allelic_freq_threshold]"
    return 1
  fi

  local input_vcf="$1"
  local min_depth="${2:-10}"  # Default minimum depth is 10 if not provided
  local allelic_freq_threshold="${3:-0.0001}"  # Default allelic frequency threshold is 0.0001 if not provided
  local output_vcf="filt_${input_vcf}"

  bcftools filter -i "INFO/DP >= ${min_depth} && INFO/POPAF < ${allelic_freq_threshold}" "$input_vcf" -o "$output_vcf"
  echo "Filtered VCF file saved as: $output_vcf"
}
