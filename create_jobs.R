#!/usr/bin/env Rscript

library(RSQLite)

main <- function() {
    script_dir <- get_script_dir()
    
    output_dir <- file.path(script_dir, 'slurm_out')
    dir.create(output_dir)
    
    jobs_dir <- file.path(script_dir, 'jobs')
    stopifnot(!file.exists(jobs_dir))
    
    jobs_file <- file(file.path(script_dir, 'jobs.txt'), 'w')
    for(use_splines in c(FALSE, TRUE)) {
        for(remove_zeros in c(FALSE, TRUE)) {
            for(use_surr_vec in list(c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))) {
                use_surr_flu <- use_surr_vec[1]
                use_surr_env <- use_surr_vec[2]
                for(use_log_flu in c(FALSE, TRUE)) {
                    for(env_var in c('T', 'AH', 'RH', 'PRCP')) {
                        for(flu_is_cause in c(FALSE, TRUE)) {
                            create_job(
                                jobs_dir,
                                jobs_file,
                                use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, env_var, flu_is_cause
                            )
                        }
                    }
                }
            }
        }
    }
    close(jobs_file)
}

create_job <- function(jobs_dir, jobs_file, use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, env_var, flu_is_cause) {
    job_dir <- file.path(jobs_dir,
        sprintf(
            'us=%d-rz=%d-usf=%d-use=%d-ulf=%d-ev=%s-fic=%d',
            use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, substr(env_var, 1, 1), flu_is_cause
        )
    )
    cat(sprintf('%s\n', job_dir), file = jobs_file)
    
    stopifnot(!file.exists(job_dir))
    dir.create(job_dir, recursive = TRUE)
    
    job_spec <- list(
        use_splines = use_splines,
        remove_zeros = remove_zeros,
        use_surr_flu = use_surr_flu,
        use_surr_env = use_surr_env,
        use_log_flu = use_log_flu,
        env_var = env_var,
        flu_is_cause = flu_is_cause
    )
    saveRDS(job_spec, file.path(job_dir, 'job_spec.Rds'))
}

get_script_dir <- function() {
    command_args <- commandArgs(trailingOnly = FALSE)
    file_arg_name <- '--file='
    tools:::file_path_as_absolute(
        dirname(
            sub(file_arg_name, '', command_args[grep(file_arg_name, command_args)])
        )
    )
}

main()
