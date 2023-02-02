# rules for running Somalier to get sex and ancestry and adding this information to the covariates file.
# Sex and first 10 PCs are used as covariates in regenie.
# # you will need an indexed vcf file for this command
rule somalier:
    """
    This takes a VCF file and extracts genotype-like information at selected sites
    Note that you can also use a BAM file as input. Please see https://github.com/brentp/somalier for more details
    @Input:
        multi/single-sample VCF from which to extract
        Sites is a VCF of known polymorphic sites in VCF format.
    @Output:
        Somalier file for each sample
    """
    input:
        vcf=input_vcf,
    output:
        somalier = expand(join(workpath, "somalier", "{sample}.somalier"), sample=SAMPLES),
        related = join(workpath, "somalier", "relatedness.samples.tsv"),
        ancestry = join(workpath, "somalier", "ancestry.somalier-ancestry.tsv"),
    params:
        rname = "somalier",
        outdir = join(workpath, "somalier"),
        ref = config['references']['GENOME'],
        sites = config['somalier']['somalier_sites'],
        ancestry_data=config['somalier']['ancestry_db'],
        exe = config['somalier']['exe'],
    shell:
        """
        echo "Extracting sites to estimate ancestry."
        {params.exe} extract -d {params.outdir}/ --sites {params.sites} -f {params.ref} {input.vcf}
        echo "Estimating relatedness."
        {params.exe} relate --infer -i -o {params.outdir}/relatedness {output.somalier}
        echo "Estimating ancestry."
        {params.exe} ancestry --n-pcs=10 -o {params.outdir}/ancestry --labels {params.ancestry_data}/ancestry-labels-1kg.tsv {params.ancestry_data}/*.somalier ++ {output.somalier} || {{
    # Somalier ancestry error,
    # usually due to not finding
    # any sites compared to the
    # its references, expected
    # with sub-sampled datasets
    echo "WARNING: Somalier ancestry failed..." 1>&2
    touch {output.ancestry}
    }}
        """

rule extract_info:
    """
    Extract sex and ancestry PCS and add to covariate file for regenie.
    @Input:
        Outputs from somalier related (5th column) and somalier ancestry (columns 9-18) tsv files
    @Output:
        Updated covariate file with sex and ancestry information
    """
    input:
        sex = join(workpath, "somalier","relatedness.samples.tsv"),
        ancestry = join(workpath, "somalier", "ancestry.somalier-ancestry.tsv"),
        covariates = covariates,
    params:
        rname = "extract_info",
        outdir_regenie=join(workpath,"regenie"),
        outdir_QC=join(workpath,"QC")
    output:
        reg_covariates = join(workpath,"regenie", "covariates.txt"),
        sex_file = join(workpath,"QC", "sex_file.txt"),
        subset = join(workpath, "somalier", "sub_ancestry_somalier.tsv")
    shell:
        """
        mkdir -p {params.outdir_regenie}
        mkdir -p {params.outdir_QC}
        paste -d '\t' <(cut -f 1,2 {input.covariates}) <(cut -f 5 {input.sex}) > {output.sex_file}
        head -n 1 {input.ancestry} > {output.subset}
        awk '/P00/' {input.ancestry} >> {output.subset}
        paste -d '\t' <(cut -f 1,2 {input.covariates}) <(cut -f 5 {input.sex}) <(cut -f 9- {output.subset}) > {output.reg_covariates}
        """
