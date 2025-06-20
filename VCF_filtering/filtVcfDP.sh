# filtVCF is a function meant to filter VCF files with a
# sequencing depth >= 10 by default
# different values can be specified in the second and third argument of the function
# the output file is saved as filt_OUTPUT

filtVcfDP() {
  if [ $#-lt 1 ] || [ $#-gt 2 ]; then
    echo "Usage: filtVcf <input_vcf_file> [min_depth]"
    return 1
  fi

  local input_vcf="$1"
  local min_depth="${2:-10}"  # Default minimum depth is 10 if not provided
  local output_vcf="filt_${input_vcf}"

  bcftools filter -i "INFO/DP >= ${min_depth}" "$input_vcf" -o "$output_vcf"
  echo "Filtered VCF file saved as: $output_vcf"
}
