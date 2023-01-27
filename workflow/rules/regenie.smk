# run step 1 and step 2 of regenie
# Result folders
regenieStep1Folder = 'working_directory/regenie/step_1'
regenieStep2Folder = 'working_directory/regenie/step_2'
regenieOutputFolder = 'working_directory/regenie/results'
docsOutputFolder = 'docs/regenie/'

rule regenie_step1:
    input:
        step1bed =
        phenoFile =
    params:
        phenoCovariates =
    conda:
        ""
    threads: 24
    shell:
        """
         regenie \
        --step 1 \
        --threads {threads} \
        --gz \
        --bed {input.step1bed} \
        """
