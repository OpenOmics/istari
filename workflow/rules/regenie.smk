# run step 1 and step 2 of regenie
# Result folders
# The chromosomes to run

rule regenie:
    input:
        step1bed = join(workpath, "regenie", 'filtered.bed'),
        phenoFile = phenotype,
        covariates = join(workpath, "regenie", 'covariates.txt'),
        step2_01 = join(workpath, "regenie", config['slivar']['gnomad_af1'] + '.bed'),
        step2_05 = join(workpath, "regenie", config['slivar']['gnomad_af2'] + '.bed'),
        annot_01 = join(workpath,"regenie", config['slivar']['gnomad_af1'] + ".annot"),
        annot_05 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".annot"),
        list_01 = join(workpath,"regenie", config['slivar']['gnomad_af1'] + ".list"),
        list_05 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".list"),
        mask_01 = join(workpath,"regenie", config['slivar']['gnomad_af1'] + ".masks"),
        mask_05 = join(workpath,"regenie", config['slivar']['gnomad_af2'] + ".masks") 
    output:
        step1 = join(workpath, "regenie", 'step1'),
        step2_01 = expand(join(workpath, "regenie", 'step2_01_{chr}.regenie'), chr=chromosomes),
        step2_05 = expand(join(workpath, "regenie", 'step2_05_{chr}.regenie'), chr=chromosomes)
    params:
        rname = "regenie",
        covarColList = config['regenie']['covarColList'],
        phenoColList = config['regenie']['phenoColList'],
        step1_bsize = config['regenie']['step1_bsize'],
        step2_bsize = config['regenie']['step2_bsize'],
    conda: config['conda']['regenie_env']
    threads: 24 
    shell:
        """
        regenie \
        --step 1 \
        --bed {input.step1bed} \
        --covarFile {input.covariates} \
        --covarColList {params.covarColList} \
        --phenoFile {input.phenoFile} \
        --phenoColList {params.phenoColList} \
        --bsize {params.step1_bsize} \
        --threads {threads} \
        --lowmem \
        --gz \
        --lowmem-prefix tmp_rg \
        --out {output.step1} \

        regenie \
        --step 2 \
        --bed {input.step2_01} \
        --covarFile {input.covariates} \
        --covarColList {params.covarColList} \
        --phenoFile {input.phenoFile} \
        --phenoColList {params.phenoColList}  \
        --write-mask-snplist \
        --pred {output.step1}_pred.list \
        --anno-file {input.annot_01} \
        --mask-def {input.mask_01} \
        --vc-tests skat,skato \
        --debug \
        --bsize {params.step2_bsize} \
        --build-mask sum \
        --check-burden-files \
        --set-list {input.list_01} \
        --strict-check-burden \
        --aaf-bins 0.01 \
        --minMAC 1 \
        --write-samples \
        --out {output.step2_01}

        regenie \
        --step 2 \
        --bed {input.step2_05} \
        --covarFile {input.covariates} \
        --covarColList {params.covarColList} \
        --phenoFile {input.phenoFile} \
        --phenoColList {params.phenoColList} \
        --write-mask-snplist \
        --pred {output.step1}_pred.list \
        --anno-file {input.annot_05} \
        --mask-def {input.mask_05} \
        --vc-tests skat,skato \
        --debug \
        --bsize {params.step2_bsize} \
        --build-mask sum \
        --check-burden-files \
        --set-list {input.list_05} \
        --strict-check-burden \
        --aaf-bins 0.05 \
        --minMAC 1 \
        --write-samples \
        --out {output.step2_05}
        """