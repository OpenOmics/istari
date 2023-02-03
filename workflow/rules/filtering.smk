# Rule for filtering for missing data using GATK, compressing, and indexing resultant file using samtools

rule GATK:
    """
    Takes merged VCF file from GLnexus Output.
    @Input:
        Multi-sample joint VCF file (scatter)
    @Output:
        VCF file filtered for missing data
    """
    input:
        vcf = input_vcf,
    output:
        cleaned_vcf = join(workpath, "QC", "cleaned.vcf.gz"),
        cleaned_vci = join(workpath, "QC", "cleaned.vcf.gz.tbi"),
        unzipped_vcf = temp(join(workpath, "QC", "cleaned.vcf")),
    params:
        rname = "GATK",
        outdir = join(workpath, "QC"),
        ref = config['references']['GENOME'],
        filt_nocall = config['QC']['GATK_filters']['max_nocall'],
    shell:
        """
        mkdir -p {params.outdir}
        set +u
        module load GATK
        gatk SelectVariants -V {input.vcf} -O {output.unzipped_vcf} -R {params.ref} --max-nocall-fraction {params.filt_nocall}
        module load bcftools
        bgzip {output.unzipped_vcf}
        tabix -p vcf {output.cleaned_vcf}
        """

rule plink:
    """
    Takes VCF file, converts to plink format and filters for missing genotype rate, removes chromosome M (regenie will produce error if included), and sets vcf half calls as missing
    @Input:
        Multi-sample joint VCF file
    @Output:
        filtered plink file for step 1 regenie, step 2 regenie files for gnomad AF choices
    """
    input:
        cleaned_vcf = join(workpath, "QC", "cleaned.vcf.gz"),
        plink_01 = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af1']),
        plink_05 = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af2'])
    params:
        rname = "plink",
        geno = config['QC']['plink_filters']['geno'],
        maf = config['QC']['plink_filters']['maf'],
        sex = join(workpath, "QC", "sex_file.txt"),
        prefix = 'geno' + config['QC']['plink_filters']['geno'] + 'maf' + config['QC']['plink_filters']['maf'],
        ld_window = config['QC']['plink_filters']['ld_window'],
        ld_step = config['QC']['plink_filters']['ld_step'],
        ld_thresh = config['QC']['plink_filters']['ld_thresh'],
        mac = config['QC']['plink_filters']['mac']
    output:
        plink1 = join(workpath, "QC",'geno' + config['QC']['plink_filters']['geno'] + 'maf' + config['QC']['plink_filters']['maf']),
        filtered = join(workpath, "regenie", 'filtered'),
        step2_01 = join(workpath, "regenie", config['slivar']['gnomad_af1']),
        step2_05 = join(workpath, "regenie", config['slivar']['gnomad_af2'])

    shell:
        """
        ml plink/1
        plink --vcf {input.cleaned_vcf} --vcf-half-call m --set-missing-var-ids @:# --not-chr M --geno {params.geno} --maf {params.maf} --update-sex {params.sex} --make-bed --out {output.plink1}
        ml plink/2
        plink2 --bfile {params.prefix} --indep-pairwise {params.ld_window} {params.ld_step} {params.ld_thresh} --make-bed --out ld_filt
        plink2 --bfile ld_filt --exclude ld_filt.prune.out --mac {params.mac} --make-bed --out {output.filtered}
        ml plink/1
        plink --bfile {output.filtered} --bmerge {input.plink_01} --make-bed --out {output.step2_01} --update-sex {params.sex}
        plink --bfile {output.filtered} --bmerge {input.plink_05} --make-bed --out {output.step2_05} --update-sex {params.sex}
        """
