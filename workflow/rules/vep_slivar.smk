# rules for annotation using VEP and slivar to subset
import os
from os.path import join
import pandas as pd

rule vep:
    """
    Takes vcf file and annotates using VEP
    @Input: multi-sample jointly-called VCF file
    @Output: VEP and HTML output
    """
    input:
        vcf = config['input']['VCF'],
    output:
        vep_vcf=join(filt_out, "cleaned.vcf.gz"),
        vep_vci=join(filt_out, "vep.vcf.gz.tbi")
    params:
        vep_assembly = config['references']['vep_assembly'],
        ref = config['references']['GENOME'],
        unzipped_vcf = temp(join(filt_out, "vep.vcf")),
   shell:
      """
      set +u
      module load VEP
      vep -i {input.vcf} -o {output.vep_vcf} --force_overwrite --fork 12 --fasta {params.ref} --species human --assembly {params.vep_assembly} --cache --dir_cache $VEP_CACHEDIR --offline --format vcf --compress_output bgzip --everything --pick --vcf
      module load bcftools
      tabix -p vcf {output.vep_vcf}
      """

rule slivar_01:
    input:
        vep_vcf = join(filt_out, "vep.vcf.gz"),
    output:
        unzipped_vcf = join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af2'] + ".vcf"),
        slivar_vcf = join(filt_out,'slivar' + 'af' + config['options']['gnomad_af1'] + ".vcf.gz"),
        slivar_annot= join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af1'] + ".txt"),
    params:
        genome = config['references']['GENOME'],
        slivar_order = config['slivar']['slivar_order'],
        slivar_tool = config['sliver']['binary'],
        slivar_js = config['sliver']['slivar_js'],
        gnomad_af = config['options']['gnomad_af1'],
   shell:
      """
      set +u
      SLIVAR_IMPACTFUL_ORDER={params.slivar_order}
      {params.slivar_tool} expr --js {params.slivar_js} -g {params.slivar_gnomad} --vcf {input.vep_vcf} --info "INFO.gnomad_popmax_af < {params.gnomad_af} && INFO.impactful" --pass-only > {output.unzipped_vcf}
      module load bcftools
      bcftools +split-vep {output.unzipped_vcf} -f '%CHROM:%POS:%REF:%ALT\t%SYMBOL\t%Consequence\n' > {output.slivar_annot}
      bgzip -c {output.unzipped_vcf} > {output.slivar_vcf}
      tabix -p vcf {output.slivar_vcf}
      """

rule slivar_05:
    input:
        vep_vcf = join(filt_out, "vep.vcf.gz"),
    output:
        unzipped_vcf = join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af2'] + ".vcf"),
        slivar_vcf = join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af2'] + ".vcf.gz"),
        slivar_vci = join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af2'] + ".vcf.gz.tbi"),
        slivar_annot= join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af2'] + ".txt"),
    params:
        genome = config['references']['GENOME'],
        slivar_order = config['slivar']['slivar_order'],
        slivar_tool = config['sliver']['binary'],
        slivar_js = config['sliver']['slivar_js'],
        gnomad_af = config['options']['gnomad_af2'],
   shell:
      """
      set +u
      SLIVAR_IMPACTFUL_ORDER={params.slivar_order}
      {params.slivar_tool} expr --js {params.slivar_js} -g {params.slivar_gnomad} --vcf {input.vep_vcf} --info "INFO.gnomad_popmax_af < {params.gnomad_af} && INFO.impactful" --pass-only > {output.unzipped_vcf}
      module load bcftools
      bcftools +split-vep {output.unzipped_vcf} -f '%CHROM:%POS:%REF:%ALT\t%SYMBOL\t%Consequence\n' > {output.slivar_annot}
      bgzip -c {output.unzipped_vcf} > {output.slivar_vcf}
      tabix -p vcf {output.slivar_vcf}
      """
rule prep_annot_files:
    """ Removes variants not within genes in annotation file
    """
    input:
        slivar_annot1= join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af1'] + ".txt"),
        slivar_annot2= join(filt_out, 'slivar' + 'af' + config['options']['gnomad_af2'] + ".txt"),
    output:
        out1 = join(workpath,"regenie", 'af' + config['options']['gnomad_af1'] + ".annot"),
        out2 = join(workpath,"regenie", 'af' + config['options']['gnomad_af2'] + ".annot")
    params:
        outdir_regenie=join(workpath,"regenie"),
 shell:
    """
    mkdir -p {params.outdir_regenie}
    awk -F'\t' '!($2 == ".")' {input.slivar_annot1} > {output.out1}
    awk -F'\t' '!($2 == ".")' {input.slivar_annot2} > {output.out2}
    dos2unix {output.out1}
    dos2unix {output.out2}
    """

rule prep_mask_list_files:
    """ Creates list file and mask file from annotation file for regenie inputs
    """
    input:
        in1 = join(workpath,"regenie", 'af' + config['options']['gnomad_af1'] + ".annot"),
        in1 = join(workpath,"regenie",'af' + config['options']['gnomad_af2'] + ".annot")
    output:
        list_out1 = join(workpath,"regenie",'af' + config['options']['gnomad_af1'] + ".list"),
        list_out2 = join(workpath,"regenie",'af' + config['options']['gnomad_af2'] + ".list"),
        mask_out1 = join(workpath,"regenie",'af' + config['options']['gnomad_af1'] + ".masks"),
        mask_out2 = join(workpath,"regenie",'af' + config['options']['gnomad_af2'] + ".masks")
    script:
        '/data/OpenOmics/istari/scipts/prepare_list_mask_files.R'
