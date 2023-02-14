# Rule for filtering for missing data using GATK, compressing, and indexing resultant file using samtools

rule GATK:
    """
    Takes merged gVCF file.
    @Input:
        Multi-sample joint VCF file (scatter)
    @Output:
        VCF file filtered for missing data
    """
    input:
        vcf = input_vcf,
    output:
        cleaned_vcf = join(workpath, "QC", "cleaned.vcf.gz"),
        cleaned_vci = join(workpath, "QC", "cleaned.vcf.gz.tbi")
    params:
        rname = "GATK",
        outdir = join(workpath, "QC"),
        unzipped_vcf = temp(join(workpath, "QC", "cleaned.vcf")),
        ref = config['references']['GENOME'],
        filt_nocall = config['QC']['GATK_filters']['max_nocall'],
    shell:
        """
        set +u
        module load GATK
        gatk SelectVariants -V {input.vcf} -O {params.unzipped_vcf} -R {params.ref} --max-nocall-fraction {params.filt_nocall}
        module load bcftools
        bgzip {params.unzipped_vcf}
        tabix -p vcf {output.cleaned_vcf}
        """

rule plink:
    """
    Takes VCF file, converts to plink format and
    filters for missing genotype rate, removes
    chromosome M (regenie will produce error if included),
    and sets vcf half calls as missing
    @Input:
        Multi-sample joint VCF file
    @Output:
        filtered plink file for step 1 regenie, step 2
        regenie files for gnomAD AF choices
    """
    input:
        cleaned_vcf = join(workpath, "QC", "cleaned.vcf.gz"),
        plink_01 = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af1']+'.bed'),
        plink_05 = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af2']+'.bed')
    params:
        rname = "plink",
        geno = config['QC']['plink_filters']['geno'],
        maf = config['QC']['plink_filters']['maf'],
        sex = join(workpath, "QC", "sex_file.txt"),
        prefix =   join(workpath, "QC",'geno' + config['QC']['plink_filters']['geno'] + 'maf' + config['QC']['plink_filters']['maf']),
        ld_window = config['QC']['plink_filters']['ld_window'],
        ld_step = config['QC']['plink_filters']['ld_step'],
        ld_thresh = config['QC']['plink_filters']['ld_thresh'],
        mac = config['QC']['plink_filters']['mac'],
        bed1 = join(workpath, "QC",'geno' + config['QC']['plink_filters']['geno'] + 'maf' + config['QC']['plink_filters']['maf']),
        ld_bed = join(workpath, "QC", "ld_filtered"),
        bed = join(workpath, "regenie", 'step1'),
        bed_01 = join(workpath, "slivar",'slivar' + config['slivar']['gnomad_af1']),
        bed_05 = join(workpath, "slivar",'slivar' + config['slivar']['gnomad_af2']),
        step2_01 = join(workpath, "regenie",'step2.' + config['slivar']['gnomad_af1']),
        step2_05 = join(workpath, "regenie", 'step2.' +config['slivar']['gnomad_af2'])
    output:
        plink1 = join(workpath, "QC",'geno' + config['QC']['plink_filters']['geno'] + 'maf' + config['QC']['plink_filters']['maf'] + '.bed'),
        ld_filt = join(workpath, "QC", "ld_filtered" + '.bed'),
        filtered = join(workpath, "regenie", 'step1' +'.bed'),
        step2_01 = join(workpath, "regenie",'step2.' + config['slivar']['gnomad_af1'] + '.bed'),
        step2_05 = join(workpath, "regenie", 'step2.' +config['slivar']['gnomad_af2'] + '.bed')

    shell:
        """
        ml plink/1
        plink --vcf {input.cleaned_vcf} --vcf-half-call m --set-missing-var-ids @:# --not-chr M --geno {params.geno} --maf {params.maf} --double-id --set-hh-missing --update-sex {params.sex} --make-bed --out {params.bed1}
        ml plink/2
        plink2 --bfile {params.prefix} --indep-pairwise {params.ld_window} {params.ld_step} {params.ld_thresh} --make-bed --out {params.ld_bed}
        plink2 --bfile {params.ld_bed} --exclude {params.ld_bed}.prune.out --mac {params.mac} --make-bed --out {params.bed}
        ml plink/1
        plink --bfile {params.bed} --bmerge {params.bed_01}  --make-bed --out {params.step2_01} --update-sex {params.sex}
        plink --bfile {params.bed} --bmerge {params.bed_05}  --make-bed --out {params.step2_05} --update-sex {params.sex}
        """
