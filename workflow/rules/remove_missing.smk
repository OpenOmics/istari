# Rule for filtering for missing data using GATK, compressing, and indexing resultant file using samtools


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
      cleaned_vcf = join(workpath, "cleaned.vcf.gz"),
      cleaned_vci = join(workpath, "cleaned.vcf.gz.tbi"),
   params:
      ref = config['references']['GENOME'],
      unzipped_vcf = temp(join(workpath, "cleaned.vcf")),
   shell:
      """
      set +u
      module load GATK/4.2.0.0
      gatk SelectVariants -V {input.vcf} -O {params.unzipped_vcf} -R {params.ref} --max-nocall-fraction 0.05
      ml samtools
      bgzip {params.unzipped_vcf}
      tabix -p vcf {output.cleaned_vcf}
      """
