# rule for running Somalier

# extract sites
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
        vcf=config['input']['VCF'],
    output:
        files="somalier_out/sites/{sample}.somalier",
    params:
        ref=config['references']['GENOME'],
        sites=config['references']['somalier_sites'],
    envmodules:
        config['tools']['somalier'],
    shell:
        "somalier extract -d {output.files} --sites {params.sites} -f {params.ref} {input.vcf}"

rule somalier_relate:
     """
    This calculates relatedness among samples from extracted, genotype-like information
    @Input:
        $sample.somalier files for each sample.
    @Output:
        This will create text and interactive HTML output that makes it fast and easy to detect mismatched samples and sample-swaps.
    """
    input:
        dir="somalier_out/sites",
    output:
        file="somalier_relate",
    envmodules:
        config['tools']['somalier'],
    shell:
        "somalier relate -i -o {output.file} {input.dir}/*.somalier"

rule somalier_ancestry:
    """
    Predict ancestry on a set of query samples
    """
    input:
        dir="somalier_out"
    output:
        file = "somalier_ancestry"
    params:
        ref=config['references']['somalier_1kg'],
    envmodules:
        config['tools']['somalier']
    shell:
        "somalier ancestry --n-pcs=10 -o {output.file} {input.dir}"

# extract first 10 pcs and create covariates file
