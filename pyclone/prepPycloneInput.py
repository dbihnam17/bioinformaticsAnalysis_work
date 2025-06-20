# Import libraries
import vcfpy
import pandas as pd

# Input VCF and CNV paths per timepoint
vcfFiles = {
    "Baseline": "/path/to/baseline_sample.muTect2.somatic.snv.vcf",
    "On_Treatment": "/path/to/on_treatment_sample.muTect2.somatic.snv.vcf",
    "Progression": "/path/to/progression_sample.muTect2.somatic.snv.vcf"
}

cnvFiles = {
    "Baseline": "/path/to/baseline_sample.final.bam_CNVs.final.gff",
    "On_Treatment": "/path/to/on_treatment_sample.final.bam_CNVs.final.gff",
    "Progression": "/path/to/progression_sample.final.bam_CNVs.final.gff"
}

sampleIDs = {
    "Baseline": "baseline_sample",
    "On_Treatment": "on_treatment_sample",
    "Progression": "progression_sample"
}

# Parse VCFs and match CNVs
rows = []

for timepoint, vcfFile in vcfFiles.items():
    sampleID = sampleIDs[timepoint]
    cnvFile = cnvFiles[timepoint]

    # Read CNV data
    cnvDf = pd.read_csv(cnvFile, sep="\t", header=None, names=[
        "chrom", "tool", "cnv_type", "start", "end", "a", "b", "c", "info"])
    cnvDf["copy_number"] = cnvDf["info"].str.extract(r"CopyNumber=(\d+)").astype(float)
    cnvDf = cnvDf[cnvDf["cnv_type"] != "breakpoint"]

    # Read VCF
    reader = vcfpy.Reader.from_path(vcfFile)
    for record in reader:
        for call in record.calls:
            if call.sample != sampleID:
                continue

            ad = call.data.get("AD")
            if not ad or len(ad) < 2:
                continue

            refCount, altCount = ad[0], ad[1]
            mutID = f"{record.CHROM}_{record.POS}_{record.REF}_{record.ALT[0].value}"

            # Match CNV for this position
            match = cnvDf[
                (cnvDf["chrom"] == record.CHROM) &
                (cnvDf["start"] <= record.POS) &
                (cnvDf["end"] >= record.POS)
            ]

            tumorCN = int(match["copy_number"].iloc[0]) if not match.empty else 2
            majorCN = tumorCN // 2 + tumorCN % 2
            minorCN = tumorCN // 2

            rows.append({
                "mutation_id": mutID,
                "sample_id": timepoint,
                "ref_counts": refCount,
                "alt_counts": altCount,
                "normal_cn": 2,
                "minor_cn": minorCN,
                "major_cn": majorCN
            })

# Save final PyClone-VI input TSV
df = pd.DataFrame(rows)

# Ensure column order
expected_columns = [
    "mutation_id", "sample_id", "ref_counts", "alt_counts",
    "normal_cn", "minor_cn", "major_cn"
    ]
df = df[expected_columns]

output_file = "sample_pycloneInput.tsv"
df.to_csv(output_file, sep="\t", index=False)
print(f"✅ PyClone-VI input saved to: {output_file}")
