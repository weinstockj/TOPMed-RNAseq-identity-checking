# Sample identity checking in TOPMed RNA-seq data

This pipeline performs two tasks:

 1. Empirically match genotypes between the RNA-seq derived genotypes and the WGS derived genotypes  
 2. Optionally run VerifyBamID on the RNA-seq files  

## Empirical genotype matching

 1. Separately by chromosome, filter the WGS derived genotypes to only include common (MAF > 5%) PASS snps. This pipeline uses the SAV export tool. 
    ```bash 
    {sav_binary} export -e 'AC>{SAV_EXPORT_MIN_AC};AC<{SAV_EXPORT_MAX_AC}' {input} |
                  bcftools view -Ob -o {output} -v snps -f PASS -

    bcftools index {output}
    ```
 2. Separately by chromosome, filter the output from 1. to the exome
    ```bash
    bcftools view -T {coding_exons} -Ob -o {output} {input}
    bcftools index {output}
    ```
 3. Concat the 22 (only autosomes) files from step 2.
    ```bash
    bcftools concat -n -Ob -o {output} {input}
    bcftools index {output}
    ```
 4. Convert output from 3. to bgzipped VCF (later tools require vcf.gz, bcf is kept for convenience)
    ```bash
    bcftools view -Oz -o {output} {input}
    bcftools index --tbi {output}
    ```
 5. Separately for each RNA-seq sample, run vt-discover2 to produce single sample BCFs
    ```bash
    {samtools_binary} view -uh -T {fasta} {params.path} 2>{log.samtools} | 
        {vt_binary} discover2 -z -q 20 -b + -r {fasta} -s {wildcards.sample_id} -o {output} 2>{log.vt}
    ```
 6. For each single sample BCF from 5., run the genotyper
    ```bash
    {vt_sample} --bcf {input} --out {output} --bed {coding_exons}
    ```
 7. For each output BCF from 6., index the bcf
    ```bash
    bcftools index {input}
    ```
 8. For each genotyped BCF from 6., empirically match the genotypes to the WGS derived common exonic vcf.gz to find the best match
    ```bash
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{htslib}
    {cramore} vcf-ibs-matrix --vcf {input.vt} --panel {input.vcf} --out {output}
    ```
 9. For each sample output in 8, print out the top 20 matches and store for convenience
    ```bash
    zcat {input} | head -n1 > {output}
    zcat {input} | tail -n+2 | sort -g -k 12 -r | head -n 20 >> {output}
    ```

## VerifyBamID

 1. Separately by chromosome, subset WGS derived genotypes to rare variants (AC <= 200)
    ```bash
    {sav_gen_info_binary} export --sites-only --generate-info SPARSE_OFFSETS_GT --filter 'AC <= 200' {input} |
    bcftools view -Ob -o {output} -f PASS -
    bcftools index {output}
    ```
 2. Subset output from 1 to exome
    ```bash
    bcftools view -T {coding_exons} -Ou {input} | bcftools sort -Ob -o {output}
    bcftools index {output}
    ```
 3. Concatenate the chromosome specific files from 2.
    ```bash
    bcftools concat -n -Ou {input} | bcftools sort -Ob -o {output}
    bcftools index {output}
    ```
 4. Use vt-discover2 to create per-sample BCF (same as step 5 in the empirical genotype matching pipeline)
 5. Identity best putative match between RNA-seq BAM files and the WGS derived rare exonic genotypes
    ```bash
    cat {sample_freeze_summary} | cut -f1 | tail -n+2 > {output.sample_ids}
    echo {input.vt} | tr ' ' '\n' | sort > bcfs_to_match.txt
    {vcf_match} --list bcfs_to_match.txt --id {output.sample_ids} --offset {input.rare_bcf} --out {output.match}
    ```
 6. Create VCF of common exonic genotypes from WGS (same as steps 1-4 in the empirical genotype matching pipepline) 
 7. Run verifyBamID on the best matching NWDID as estimated from variants using WGS derived common exonic genotypes as the reference VCF
    ```bash
    best_match=$(tail -n+2 {input.match} | awk '{{if ($1 == "{wildcards.sample_id}") {{print $3}} }}' | cut -d":" -f1)
    echo "best match for {wildcards.sample_id} is $best_match"
    {verifyBamID} --smID $best_match --vcf {input.vcf} --bam {params.path} --out {verify_dir}/{wildcards.sample_id} --self --ignoreRG
    ```

## Software details

Here is a list of CSG specific hardcoded paths to the required binaries (making this more portable would be an improvement):
```bash
sav_gen_info_binary = "/net/wonderland/home/lefaivej/savvy-dev/build/sav-gen-info"
sav_binary = "/net/fantasia/home/hmkang/tools/savvy/bin/sav"

samtools_binary = "/net/topmed8/working/call_sets/freeze9/topmed_variant_calling/samtools/samtools"
vt_binary = "/net/topmed8/working/call_sets/freeze9/topmed_variant_calling/vt-topmed/vt"

coding_exons = "/net/1000g/hmkang/data/gencode/gtfs/gencode.v34.annotation.gtf.coding_exons.sorted_merged.bed"

fasta = "/net/topmed8/working/call_sets/rnaseqQC/ref/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"

vcf_match = "/net/fantasia/home/hmkang/tools/apigenome.master/bin/vcf-match-sparse-offsets"
vt_sample = "/net/fantasia/home/hmkang/tools/apigenome.master/bin/vcf-gt-vt-sample"

verifyBamID = "/net/fantasia/home/hmkang/code/working/verifyBamID/bin/verifyBamID"

cramore = "/net/fantasia/home/hmkang/code/working/cramore/bin/cramore"
htslib = "/net/fantasia/home/hmkang/code/working/htslib"
```
### Snakemake details
Snakemake 5.10.0 and pandas 0.25.1 are required. Creating a conda environment would represent an improvement here. 
I recommend installation with pip or conda. I run snakemake with the following command:

```bash
set -eou pipefail

num_jobs=300
snakemake --cluster-config cluster.yaml --cluster-sync "srun --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus} --partition {cluster.partition} --output {cluster.output} --error {cluster.error} --job-name {cluster.name}" -j $num_jobs --keep-going
```

## Authors

Josh Weinstock (jweinstk@umich.edu) and Hyun Min Kang (hmkang@umich.edu). Hyun contributed most of the 
underlying scripts and command line tools, and Josh wrapped them into a pipeline. 
