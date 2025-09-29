rule log_session_info:
    """
    Log R session information and package versions from the analysis environment
    """
    input:
        # Wait for key analysis outputs to ensure packages are loaded
        expand(OUTPUT_CANCER, cancer=CANCER_TYPES),
        BATCH_CORRECTED_EXPRESSION_FILE,
        PORCUPINE_RESULTS_ALL
    output:
        r_session = os.path.join(OUTPUT_DIR_ALL_CANCERS, "logs", "R_session_info.txt")
    shell:
        """
        mkdir -p $(dirname {output.r_session})
        
        Rscript -e "
        cat('=== R SESSION INFORMATION ===\\n', file='{output.r_session}')
        cat('Date:', as.character(Sys.time()), '\\n', file='{output.r_session}', append=TRUE)
        cat('R Version:', R.version.string, '\\n', file='{output.r_session}', append=TRUE)
        cat('Platform:', R.version\$platform, '\\n', file='{output.r_session}', append=TRUE)
        cat('\\n=== SESSION INFO ===\\n', file='{output.r_session}', append=TRUE)
        capture.output(sessionInfo(), file='{output.r_session}', append=TRUE)
        cat('\\n=== INSTALLED PACKAGES ===\\n', file='{output.r_session}', append=TRUE)
        installed <- installed.packages()
        pkg_df <- data.frame(Package=installed[,'Package'], Version=installed[,'Version'])
        write.table(pkg_df, '{output.r_session}', sep='\\t', row.names=FALSE, quote=FALSE, append=TRUE)
        "
        """