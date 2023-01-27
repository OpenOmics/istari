# rules for running Somalier to get sex and ancestry and adding this information to the covariates file.
# Sex and first 10 PCs are used as covariates in regenie.

rule somalier_extract:
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
        files=expand(join(workpath, "somalier_output", "sites", "{sample}.somalier"), sample=SAMPLES)
    params:
        ref=config['references']['GENOME'],
        sites=config['somalier']['somalier_sites'],
        exe=config['somalier']['som_software'],
    shell:
        """
        {params.exe} extract -d {output.files} --sites {params.sites} -f {params.ref} {input.vcf}
        """


rule somalier_relate:
    """
    Calculate relatedness on the extracted data. This will give you sex information for individuals.
    @Input:
        Somalier *sites files from somalier_extract
    @Output:
        Text and interactive HTML output.
    """
    input:
        sites = expand(join(workpath, "somalier_output", "sites", "{sample}.somalier"), sample=SAMPLES)
    output:
        file = join(workpath, "somalier_output", "somalier_relate.samples.tsv")
    params:
        exe=config['somalier']['som_software'],
        prefix=join(workpath, "somalier_output", "somalier_relate")
    shell:
        """
        {params.exe} relate -i -o {params.prefix} {input.sites}
        """


rule somalier_ancestry:
    """
    Perform ancestry prediction on a set of samples, given a set of labeled samples
        @Input:
            A set of labelled samples
        @Output:
            This command will create an html output along with a text file of the predictions.
    """
    input:
        dir= expand(join(workpath, "somalier_output", "sites", "{sample}.somalier"), sample=SAMPLES)
    output:
        file = join(workpath, "somalier_output", "somalier_ancestry.somalier-ancestry.tsv")
    params:
        ref=config['somalier']['somalier_1kg'],
        exe=config['somalier']['som_software'],
        prefix=join(workpath, "somalier_output", "somalier_ancestry")
    shell:
        """
        {params.exe} ancestry --n-pcs=10 -o {params.prefix} --labels {params.ref} ++ {input.dir}
        """

rule covar_file:
    """
    Extract sex and ancestry PCS and add to covariate file for regenie.
    @Input:
        Outputs from somalier related (5th column) and somalier ancestry (columns 9-18) tsv files
    @Output:
        Updated covariate file with sex and ancestry information
    """
    input:
        sex = join(workpath, "somalier_output","somalier_relate.samples.tsv"),
        ancestry = join(workpath, "somalier_output", "somalier_ancestry.somalier-ancestry.tsv"),
        covariates = covariates,
    params:
        outdir_regenie=join(workpath,"regenie")
    output:
        reg_covariates = join(workpath,"regenie", "covariates.txt")
    shell:
        """
        mkdir -p {params.outdir_regenie}
        dos2unix {input.covariates}
        dos2unix {input.sex}
        dos2unix {input.ancestry}
        paste -d '\t' <(cut -f 1,2 {input.covariates}) <(cut -f 5 {input.sex}) <(cut -f 9- {input.ancestry}) > {output.reg_covariates}
        """
