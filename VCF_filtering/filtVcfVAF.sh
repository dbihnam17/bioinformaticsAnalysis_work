# filtVCF is a function meant to filter VCF files with a
# VAF above a specified number
# The output file is stored as VAFfilt_{min_vaf}_{input_vcf}

filtVcfVAF() {
  if [ $#-lt 1 ] || [ $#-gt 2 ]; then
    echo "Usage: filtVcf <input_vcf_file> [min_VAF]"
    return 1
  fi

  if ! command -v bcftools &> /dev/null; then
    echo "bcftools not found, attempting to load module..."
    module load bcftools || { echo "Failed to load bcftools"; return 1; }
  fi

  local input_vcf="$1"
  local min_vaf="${2:-.10}"  # Default VAF is .01 if not provided
  local min_vaf_int=$(echo "$min_vaf" | awk '{printf "%d", $1 * 100}')
  local output_vcf="VAFfilt_${min_vaf_int}_${input_vcf}"

  bcftools filter -i "FORMAT/AF >= ${min_vaf}" "$input_vcf" -o "$output_vcf"
  echo "Filtered VCF file saved as: $output_vcf"
}
