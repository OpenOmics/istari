# Rule for filtering for missing data using GATK, compressing, and indexing resultant file using samtools

# fix workpath
rule remove_missing:
     """
    Takes merged VCF file from GLnexus Output.
    @Input:
        Multi-sample joint VCF file (scatter)
    @Output:
        VCF file filtered for missing data
    """
    input:
        vcf = config['input']['VCF'],
    output:
        cleaned_vcf = join(filt_out, "cleaned.vcf.gz"),
        cleaned_vci = join(filt_out, "cleaned.vcf.gz.tbi"),
    params:

        ref = config['references']['GENOME'],
        filt_nocall = config['QC']['max_nocall'],
        unzipped_vcf = temp(join(filt_out, "cleaned.vcf")),
   shell:
      """
      mkdir -p {params.outdir}
      set +u
      module load GATK/4.2.0.0
      gatk SelectVariants -V {input.vcf} -O {params.unzipped_vcf} -R {params.ref} --max-nocall-fraction {params.filt_nocall}
      ml samtools
      bgzip {params.unzipped_vcf}
      tabix -p vcf {output.cleaned_vcf}
      """
# something with thinning snps to 999999?
rule plink_convert:
     """
    Takes VCF file, converts to plink format and filters for missing genotype rate, removes chromosome M (regenie will produce error if included), and sets vcf half calls as missing
    @Input:
        Multi-sample joint VCF file
    @Output:
        filtered plink file
    """
    input:
        vcf = cleaned_vcf = join(filt_out, "cleaned.vcf.gz"),
    params:
        filt_geno = config['QC']['geno'],
    output:
        plink1 = join(filt_out, "geno_filt",config['QC']['geno'] + '.bed')
    shell:
    """
    ml plink/1
    plink --vcf {input.vcf} --vcf-half-call m --set-missing-var-ids @:# --not-chr M --geno 0.2 --make-bed --out {output.plink1}
    """
