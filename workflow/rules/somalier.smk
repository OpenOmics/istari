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
        vcf=config['input']['VCF']
    output:
        files=join(somalier_output, "sites", "{sample}.somalier")
    params:
        ref=config['references']['GENOME'],
        sites=config['somalier']['somalier_sites']
    envmodules:
        config['somalier']['som_software']
    shell:
    """
    somalier extract -d {output.files} --sites {params.sites} -f {params.ref} {input.vcf}
    """

rule somalier_relate
    """
    Calculate relatedness on the extracted data. This will give you sex information for individuals.
    @Input:
        Somalier *sites files from somalier_extract
    @Output:
        Text and interactive HTML output.
    """
    input:
        sites = join(somalier_output, "sites", "{sample}.somalier")
    envmodules:
        config['somalier']['som_software']
    output:
        file = join(somalier_output,"somalier_relate")
    shell:
    """
    somalier relate -i -o {output.file} {input.sites}
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
        dir= join(somalier_output, "sites", "{sample}.somalier")
    output:
        file = join(somalier_output,"somalier_ancestry")
    params:
        ref=config['somalier']['somalier_1kg']
    envmodules:
        config['somalier']['som_software']
    shell:
    """
    somalier ancestry --n-pcs=10 -o {output.file} --labels {params.ref} ++ {input.dir}
    """

rule covar_file
    """
    Extract sex and ancestry PCS and add to covariate file for regenie.
    @Input:
        Outputs from somalier related (5th column) and somalier ancestry (columns 9-18) tsv files
    @Output:
        Updated covariate file with sex and ancestry information
    """
    input:
        sex = join(somalier_out,"somalier_relate.samples.tsv"),
        ancestry = join(somalier_out, "/somalier_ancestry.somalier-ancestry.tsv"),
        covariates = config['input']['covarFile']
    shell:
    """
        dos2unix {input.covariates}
        dos2unix {input.sex}
        dos2unix {input.ancestry}
        paste -d '\t' <(cut -f 1,2 {input.covariates}) <(cut -f 5 {input.sex}) <(cut -f 9- {input.ancestry}) > regenie_covariates.txt
    """
