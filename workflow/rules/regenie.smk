# run step 1 and step 2 of regenie
# Result folders
# Run each chromosome separately (recommended for regenie step2 with many variants)

rule regenie:
    input:
        annot_01 = join(workpath,"regenie", 'slivar' + config['slivar']['gnomad_af1'] + ".annot"),
        annot_05 = join(workpath,"regenie", 'slivar' + config['slivar']['gnomad_af2'] + ".annot"),
        list_01 = join(workpath,"regenie", 'slivar' + config['slivar']['gnomad_af1'] + ".list"),
        list_05 = join(workpath,"regenie", 'slivar' + config['slivar']['gnomad_af2'] + ".list"),
        mask_01 = join(workpath,"regenie", 'slivar' + config['slivar']['gnomad_af1'] + ".masks"),
        mask_05 = join(workpath,"regenie", 'slivar' + config['slivar']['gnomad_af2'] + ".masks")
    output:
        step1out = join(workpath,"regenie", 'step1_pred.list'),
        step2_01 = expand(join(workpath, "regenie", 'step2_01_{pheno}.regenie'), pheno=phenos),
        step2_05 = expand(join(workpath, "regenie", 'step2_05_{pheno}.regenie'), pheno=phenos)
    params:
        rname = "regenie",
        step1 = join(workpath, "regenie", 'step1'),
        covariates = join(workpath, "regenie", 'covariates.txt'),
        covarColList = config['regenie']['covarColList'],
        phenoFile = phenotype,
        phenoColList = config['regenie']['phenoColList'],
        step1_bsize = config['regenie']['step1_bsize'],
        bed_01 = join(workpath, "regenie",'step2.'+ config['slivar']['gnomad_af1']),
        bed_05 = join(workpath, "regenie", 'step2.'+ config['slivar']['gnomad_af2']),
        step2_bsize = config['regenie']['step2_bsize'],
        step201 = join(workpath, "regenie", 'step2_01'),
        step205 = join(workpath, "regenie", 'step2_05')
    conda: config['conda']['regenie_env']
    threads: 24
    shell:
        """
        regenie \
        --step 1 \
        --bed {params.step1} \
        --covarFile {params.covariates} \
        --covarColList {params.covarColList} \
        --phenoFile {params.phenoFile} \
        --phenoColList {params.phenoColList} \
        --bsize {params.step1_bsize} \
        --threads {threads} \
        --lowmem \
        --gz \
        --nauto 23 \
        --lowmem-prefix tmp_rg \
        --out {params.step1} \

        regenie \
        --step 2 \
        --bed {params.bed_01} \
        --covarFile {params.covariates} \
        --covarColList {params.covarColList} \
        --phenoFile {params.phenoFile} \
        --phenoColList {params.phenoColList} \
        --write-mask-snplist \
        --pred {params.step1}_pred.list \
        --anno-file {input.annot_01} \
        --mask-def {input.mask_01} \
        --vc-tests skat,skato \
        --debug \
        --nauto 23 \
        --bsize {params.step2_bsize} \
        --build-mask sum \
        --check-burden-files \
        --set-list {input.list_01} \
        --strict-check-burden \
        --aaf-bins 0.01 \
        --minMAC 1 \
        --write-samples \
        --out {params.step201} \

        regenie \
        --step 2 \
        --bed {params.bed_05} \
        --covarFile {params.covariates} \
        --covarColList {params.covarColList} \
        --phenoFile {params.phenoFile} \
        --phenoColList {params.phenoColList} \
        --write-mask-snplist \
        --pred {params.step1}_pred.list \
        --anno-file {input.annot_05} \
        --mask-def {input.mask_05} \
        --vc-tests skat,skato \
        --debug \
        --nauto 23 \
        --bsize {params.step2_bsize} \
        --build-mask sum \
        --check-burden-files \
        --set-list {input.list_05} \
        --strict-check-burden \
        --aaf-bins 0.05 \
        --minMAC 1 \
        --write-samples \
        --out {params.step205} \
        """
