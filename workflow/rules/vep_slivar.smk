# rules for annotation using VEP and slivar to subset and creating mask, list, and annotation files for gene-based tests in regenie

rule vep:
    """
    Takes vcf file, filter out anything
    but biallelic SNPs and annotates using VEP
    @Input: multi-sample jointly-called VCF file
    @Output: VEP and HTML output
    """
    input:
        vcf = input_vcf
    output:
        sub = join(workpath, "vep", "2alleles.vcf.gz")
        vep_vcf=join(workpath,"vep", "vep.vcf.gz"),
        vep_vci=join(workpath,"vep", "vep.vcf.gz.tbi")
    params:
        rname = "vep",
        outdir = (workpath,"vep"),
        vep_assembly = config['references']['vep_assembly'],
        ref = config['references']['GENOME'],
    shell:
        """
        module load bcftools
        bcftools view --max-alleles 2 {input.vcf} > {output.sub}
        mkdir -p {params.outdir}
        set +u
        module load VEP
        vep -i {output.sub} -o {output.vep_vcf} --force_overwrite --fork 12 --fasta {params.ref} --species human --assembly {params.vep_assembly} --cache --dir_cache $VEP_CACHEDIR --offline --format vcf --compress_output bgzip --everything --pick --vcf
        module load bcftools
        tabix -p vcf {output.vep_vcf}
        """

rule slivar_01:
    input:
        vep_vcf = join(workpath, "vep", "vep.vcf.gz")
    output:
        unzipped_vcf = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af1'] + ".vcf"),
        slivar_vcf = join(workpath, "slivar",'slivar' + config['slivar']['gnomad_af1'] + ".vcf.gz"),
        slivar_vci = join(workpath, "slivar", 'slivar'  + config['slivar']['gnomad_af1'] + ".vcf.gz.tbi"),
        slivar_annot= join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af1'] + ".txt"),
        plink_out = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af1'])
    params:
        rname = "slivar_01",
        outdir = join(workpath,"slivar"),
        gnomad_db = config['slivar']['gnomad_db'],
        genome = config['references']['GENOME'],
        slivar_order = config['slivar']['slivar_order'],
        slivar_tool = config['slivar']['binary'],
        slivar_js = config['slivar']['slivar_js'],
        gnomad_af = config['slivar']['gnomad_af1']
    shell:
        """
        mkdir -p {params.outdir}
        set +u
        SLIVAR_IMPACTFUL_ORDER={params.slivar_order}
        {params.slivar_tool} expr --js {params.slivar_js} -g {params.gnomad_db} --vcf {input.vep_vcf} --info "INFO.gnomad_popmax_af < {params.gnomad_af} && INFO.impactful" --pass-only > {output.unzipped_vcf}
        module load plink/2
        plink2 --make-bed --max-alleles 2 --out {output.plink_out} --rm-dup force-first --set-missing-var-ids @:# --vcf {output.unzipped_vcf} --vcf-half-call m
        module load bcftools
        bcftools +split-vep {output.unzipped_vcf} -f '%CHROM:%POS:%REF:%ALT\t%SYMBOL\t%Consequence\n' > {output.slivar_annot}
        bgzip -c {output.unzipped_vcf} > {output.slivar_vcf}
        tabix -p vcf {output.slivar_vcf}
        """

rule slivar_05:
    input:
        vep_vcf = join(workpath, "vep", "vep.vcf.gz")
    output:
        unzipped_vcf = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af2'] + ".vcf"),
        slivar_vcf = join(workpath, "slivar",'slivar' + config['slivar']['gnomad_af2'] + ".vcf.gz"),
        slivar_vci = join(workpath, "slivar", 'slivar'  + config['slivar']['gnomad_af2'] + ".vcf.gz.tbi"),
        slivar_annot= join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af2'] + ".txt"),
        plink_out = join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af2'])
    params:
        rname = "slivar_05",
        outdir = join(workpath,"slivar"),
        gnomad_db = config['slivar']['gnomad_db'],
        genome = config['references']['GENOME'],
        slivar_order = config['slivar']['slivar_order'],
        slivar_tool = config['slivar']['binary'],
        slivar_js = config['slivar']['slivar_js'],
        gnomad_af = config['slivar']['gnomad_af2']
    shell:
        """
        mkdir -p {params.outdir}
        set +u
        SLIVAR_IMPACTFUL_ORDER={params.slivar_order}
        {params.slivar_tool} expr --js {params.slivar_js} -g {params.gnomad_db} --vcf {input.vep_vcf} --info "INFO.gnomad_popmax_af < {params.gnomad_af} && INFO.impactful" --pass-only > {output.unzipped_vcf}
        module load plink/2
        plink2 --make-bed --max-alleles 2 --out {output.plink_out} --rm-dup force-first --set-missing-var-ids @:# --vcf {output.unzipped_vcf} --vcf-half-call m
        module load bcftools
        bcftools +split-vep {output.unzipped_vcf} -f '%CHROM:%POS:%REF:%ALT\t%SYMBOL\t%Consequence\n' > {output.slivar_annot}
        bgzip -c {output.unzipped_vcf} > {output.slivar_vcf}
        tabix -p vcf {output.slivar_vcf}
        """

rule prep_annot_files:
    """ Removes variants not within genes in annotation file
    """
    input:
        slivar_annot1= join(workpath, "slivar", 'slivar' + config['slivar']['gnomad_af1'] + ".txt"),
        slivar_annot2= join(workpath, "slivar", 'slivar'  + config['slivar']['gnomad_af2'] + ".txt")
    output:
        out1 = join(workpath,"regenie", config['slivar']['gnomad_af1'] + ".annot"),
        out2 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".annot")
    params:
        rname = "prep_annot_files",
    shell:
        """
        awk -F'\t' '!($2 == ".")' {input.slivar_annot1} > {output.out1}
        awk -F'\t' '!($2 == ".")' {input.slivar_annot2} > {output.out2}
        dos2unix {output.out1}
        dos2unix {output.out2}
        """

rule regenie_files:
    """ Creates list file and mask file from annotation file for regenie inputs
    """
    input:
        in1 = join(workpath,"regenie",  config['slivar']['gnomad_af1'] + ".annot"),
        in2 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".annot")
    output:
        list_out1 = join(workpath,"regenie", config['slivar']['gnomad_af1'] + ".list"),
        list_out2 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".list"),
        mask_out1 = join(workpath,"regenie", config['slivar']['gnomad_af1'] + ".masks"),
        mask_out2 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".masks")
    script:
        '/data/OpenOmics/istari/scipts/prepare_list_mask_files.R'
