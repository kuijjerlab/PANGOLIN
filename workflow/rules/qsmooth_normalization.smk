# PySNAIL qsmooth normalization of combined expression data
rule qsmooth_normalization:
    """
    Apply qsmooth normalization to combined TCGA expression data using PySNAIL.
    
    Prerequisites: PySNAIL package must be installed in envs/ directory.
    Installation: git clone git@github.com:kuijjerlab/PySNAIL.git

    """
    input:
        xprs = OUTPUT_EXP_COMBINED_FILE,
        groups = GROUP_FILE
    output:
        norm = PYSNAIL_NORMALIZED_FILE
    log:
        "logs/qsmooth_normalization.log"
    message:
        "Performing qsmooth normalization on combined expression data"
    params:
        threshold = 0.2,
        bin = config["bin"]
    conda:
        PYSNAIL_YAML_RELATIVE 
    shell:
        """
        set +u
        unset PYTHONPATH
        unset PYTHONHOME
        export PYTHONNOUSERSITE=1
        python {params.bin}/normalize_with_pysnail.py {input.xprs} {input.groups} {output.norm} --threshold {params.threshold} \
        > {log} 2>&1
        """

rule split_expression_by_cancer:
    """
    Split the large normalized expression file into cancer-specific files.
    """
    input:
        expression_file = PYSNAIL_NORMALIZED_FILE,
        group_file = GROUP_FILE
    output:
        pysnail_normalized_file_cancer_specific = PYSNAIL_NORMALIZED_FILE_CANCER_SPECIFIC
    log:
        "logs/split_expression_by_{cancer}.log"
    message:
        "Splitting expression data into cancer-specific files and saving them: {wildcards.cancer}"
    params:
        bin = config["bin"]
    shell:
        """
        Rscript {params.bin}/save_normalized_exp_per_cancer.R \
            --expression_file {input.expression_file} \
            --group_file {input.group_file} \
            --cancer {wildcards.cancer} \
            --output_file {output.pysnail_normalized_file_cancer_specific} \
        > {log} 2>&1
        """

