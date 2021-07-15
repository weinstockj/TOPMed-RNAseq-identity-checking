################################################################################
# Author: Josh Weinstock <jweinstk@umich.edu>, Hyun Min Kang <hmkang@umich.edu>
# Sample identity checking in TOPMed RNA-seq
# Last updated 7/15/2021
################################################################################

import os
import pandas as pd

########### Change these parameters to alter input sample list, output name, or name of freeze #######
# expects two columns (without headers) - the first with sample name, second absolute path
# manifest_file = "/net/topmed8/working/call_sets/rnaseqQC/metadata/broad_uw_bamlist.rand100.index"
manifest_file = "/net/topmed8/working/call_sets/rnaseqQC/metadata/broad_uw_bamlist.index"
match_output_filename = "match.rand100.out"
freeze_label = "9"
sample_freeze_summary=f"/net/topmed8/working/call_sets/freeze{freeze_label}/release/summary/topmed/samples/freeze.{freeze_label}.chr22.pass_and_fail.gtonly.minDP0.sample_summary"
freeze = f"freeze{freeze_label}"

# whi only samples for debugging
# whi_file = "/net/topmed2/working/jweinstk/identity_checking_rnaseq/nhlbi.6607.WHI.rna.seq.samples.TOR"

SAV_EXPORT_MIN_AC=16097 # MAF of 5% when N = 160,970
SAV_EXPORT_MAX_AC=305851 # MAF of 95% when N = 160,970 
chroms = list(map(lambda x: f"chr{x}", range(1, 23))) # run only on autosomes

#####################################################################################################

########## Binaries and data files required by the pipeline ########
sav_dir = f"/net/topmed8/working/call_sets/freeze{freeze_label}/release/summary/topmed/sav/"
# chroms, = glob_wildcards(os.path.join(sav_dir, "freeze.9.{id}.pass_and_fail.gtonly.minDP0.sav"))

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
####################################################################################################

sample_manifest = pd.read_table(manifest_file, header = None, names = ['sample_id', 'path'])
sample_ids = sample_manifest.sample_id.values
paths = sample_manifest.path.values
# print(f"running on {len(paths)} samples")

# whi_sample_ids = []
# with open(whi_file, "r") as f:
#     for line in f:
#         whi_sample_ids.append(line.strip())

###### Output directory names ################################
rare_subset_dir = os.path.join(freeze, "rare_variant_subset")
rare_subset_exons_dir = os.path.join(freeze, "rare_variant_subset_exons")
common_subset_dir = os.path.join(freeze, "common_variant_subset")
common_subset_exons_dir = os.path.join(freeze, "common_variant_subset_exons")
vt_dir = os.path.join(freeze, "vt")
vt_gt_dir = os.path.join(freeze, "vt_gt")
cohort_vt_gt_dir = os.path.join(freeze, "cohort_vt_gt")
log_dir = os.path.join(freeze, "logs")
match_dir = os.path.join(freeze, "match")
ibs_matrix_dir = os.path.join(freeze, "ibs_matrix")
ibs_matrix_head_dir = os.path.join(freeze, "ibs_matrix_head")
verify_dir = os.path.join(freeze, "verifyBamID")
#######################################################

rule all:
    input:
        os.path.join(rare_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.maxac200.offsets.exons.bcf"),
        os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.minmaf05.exons.bcf"),
        expand(os.path.join(vt_dir, "{sample_id}.bcf"), sample_id = sample_ids),
        expand(os.path.join(vt_gt_dir, "{sample_id}.gt.bcf"), sample_id = sample_ids),
        expand(os.path.join(vt_gt_dir, "{sample_id}.gt.bcf.csi"), sample_id = sample_ids),
        expand(os.path.join(ibs_matrix_dir, "{sample_id}.txt.gz"), sample_id = sample_ids),
        expand(os.path.join(ibs_matrix_head_dir, "{sample_id}.head.txt"), sample_id = sample_ids)
        # os.path.join(cohort_vt_gt_dir, "whi.rnaseq.gt.bcf")
        # os.path.join(match_dir, match_output_filename),
        # expand(os.path.join(verify_dir, "{sample_id}.selfSM"), sample_id = sample_ids)

rule subset_to_rare:
    input:
        os.path.join(sav_dir, f"freeze.{freeze_label}.{{chrom}}.pass_and_fail.gtonly.minDP0.sav"),
    output:
        os.path.join(rare_subset_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.maxac200.offsets.bcf")
    shell:
        """
            {sav_gen_info_binary} export --sites-only --generate-info SPARSE_OFFSETS_GT --filter 'AC <= 200' {input} |
            bcftools view -Ob -o {output} -f PASS -
            bcftools index {output}
        """
        
rule subset_rare_variants_to_exons:
    input: 
        os.path.join(rare_subset_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.maxac200.offsets.bcf")
    output:
        os.path.join(rare_subset_exons_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.maxac200.offsets.exons.bcf")
    shell:
        """
            bcftools view -T {coding_exons} -Ou {input} | bcftools sort -Ob -o {output}
            bcftools index {output}
        """

rule concat_rare_coding_variants:
    input:
        expand(os.path.join(rare_subset_exons_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.maxac200.offsets.exons.bcf"), chrom = chroms)
    output:
        os.path.join(rare_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.maxac200.offsets.exons.bcf")
    shell:
        """
            bcftools concat -n -Ou {input} | bcftools sort -Ob -o {output}
            bcftools index {output}
        """

rule subset_common:
    input:
        os.path.join(sav_dir, f"freeze.{freeze_label}.{{chrom}}.pass_and_fail.gtonly.minDP0.sav")
    output:
        os.path.join(common_subset_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.minmaf05.bcf")
    shell:
        """
            # {sav_binary} export -e 'AC>16097;AC<305851' {input} | 
            {sav_binary} export -e 'AC>{SAV_EXPORT_MIN_AC};AC<{SAV_EXPORT_MAX_AC}' {input} | 
            bcftools view -Ob -o {output} -v snps -f PASS -
            bcftools index {output}
        """

rule subset_common_variants_to_exons:
    input:
        os.path.join(common_subset_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.minmaf05.bcf")
    output:
        os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.minmaf05.exons.bcf")
    shell:
        """
            bcftools view -T {coding_exons} -Ob -o {output} {input}
            bcftools index {output}
        """
    
rule concat_common_coding_variants:
    input:
        expand(os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.{{chrom}}.pass.gtonly.minDP0.minmaf05.exons.bcf"), chrom = chroms)
    output:
        os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.minmaf05.exons.bcf")
    shell:
        """
            bcftools concat -n -Ob -o {output} {input}
            bcftools index {output}
        """
rule convert_common_coding:
    input:
        os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.minmaf05.exons.bcf")
    output:
        os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.minmaf05.exons.vcf.gz")
    shell:
        """
        bcftools view -Oz -o {output} {input}
        bcftools index --tbi {output}
        """

rule vt:
    output:
        os.path.join(vt_dir, "{sample_id}.bcf")
    params:
        path = lambda wildcards: sample_manifest.path[sample_manifest.sample_id == wildcards.sample_id].values[0]
    threads: 2
    log: 
        samtools = os.path.join(log_dir, "samtools_{sample_id}.stderr.log"),
        vt = os.path.join(log_dir, "vt_{sample_id}.stderr.log")
    shell:
        """
        {samtools_binary} view -uh -T {fasta} {params.path} 2>{log.samtools} | {vt_binary} discover2 -z -q 20 -b + -r {fasta} -s {wildcards.sample_id} -o {output} 2>{log.vt}
        """
        
rule best_match_with_rare_variants:
    input:
        rare_bcf = os.path.join(rare_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.maxac200.offsets.exons.bcf"),
        vt = expand(os.path.join(vt_dir, "{sample_id}.bcf"), sample_id = sample_ids)
    output:
        sample_ids = os.path.join(rare_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.maxac200.offsets.exons.bcf.sample_ids"),
        match = os.path.join(match_dir, match_output_filename)
    shell:
        """
        cat {sample_freeze_summary} | cut -f1 | tail -n+2 > {output.sample_ids}
        echo {input.vt} | tr ' ' '\n' | sort > bcfs_to_match.txt
        {vcf_match} --list bcfs_to_match.txt --id {output.sample_ids} --offset {input.rare_bcf} --out {output.match}
        """

rule vt_genotype:
    input:
        os.path.join(vt_dir, "{sample_id}.bcf")
    output:
        vt = os.path.join(vt_gt_dir, "{sample_id}.gt.bcf")
    shell:
        """
        {vt_sample} --bcf {input} --out {output} --bed {coding_exons}
        """

rule vt_genotype_index:
    input:
        os.path.join(vt_gt_dir, "{sample_id}.gt.bcf")
    output:
        os.path.join(vt_gt_dir, "{sample_id}.gt.bcf.csi")
    shell:
        """
            bcftools index {input}
        """

# rule merge_vt_genotype_whi:
#     input:
#         bcf = expand(os.path.join(vt_gt_dir, "{sample_id}.gt.bcf"), sample_id = whi_sample_ids),
#         csi = expand(os.path.join(vt_gt_dir, "{sample_id}.gt.bcf.csi"), sample_id = sample_ids)
#     output:
#         os.path.join(cohort_vt_gt_dir, "whi.rnaseq.gt.bcf")
#     shell:
#         """
#             bcftools merge -Ob -o {output} {input.bcf}
#             bcftools index {output}
#         """

rule best_match_with_common_variants:
    input:
        vcf = os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.minmaf05.exons.vcf.gz"),
        vt = os.path.join(vt_gt_dir, "{sample_id}.gt.bcf")
    output:
        os.path.join(ibs_matrix_dir, "{sample_id}.txt.gz")
    shell:
        """
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{htslib}
        {cramore} vcf-ibs-matrix --vcf {input.vt} --panel {input.vcf} --out {output}
        """

rule ibs_matrix_head:
    input:
        os.path.join(ibs_matrix_dir, "{sample_id}.txt.gz")
    output:
        os.path.join(ibs_matrix_head_dir, "{sample_id}.head.txt")
    shell:
        """
        set +o pipefail
        set +u
        zcat {input} | head -n1 > {output}
        zcat {input} | tail -n+2 | sort -g -k 12 -r | head -n 20 >> {output}
        """
        
rule verifyBamID:
    input:
        match = os.path.join(match_dir, match_output_filename),
        vcf = os.path.join(common_subset_exons_dir, f"freeze.{freeze_label}.autosomes.pass.gtonly.minDP0.minmaf05.exons.vcf.gz")
    params:
        path = lambda wildcards: sample_manifest.path[sample_manifest.sample_id == wildcards.sample_id].values[0]
    output:
        os.path.join(verify_dir, "{sample_id}.selfSM")
    shell:
        """
        best_match=$(tail -n+2 {input.match} | awk '{{if ($1 == "{wildcards.sample_id}") {{print $3}} }}' | cut -d":" -f1)
        echo "best match for {wildcards.sample_id} is $best_match"
        {verifyBamID} --smID $best_match --vcf {input.vcf} --bam {params.path} --out {verify_dir}/{wildcards.sample_id} --self --ignoreRG
        """
