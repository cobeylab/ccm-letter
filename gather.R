#!/usr/bin/env Rscript

library(RSQLite)

main <- function() {
    script_dir <- get_script_dir()
    
    outdb_filename <- file.path(script_dir, 'output_db_all.sqlite')
    stopifnot(!file.exists(outdb_filename))
    
    db <- dbConnect(SQLite(), outdb_filename)
    dbGetQuery(db, '
        CREATE TABLE ccm_rho (
            use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, cause, effect, country,
            rho_sig_95, rho, rho_null_025, rho_null_05, rho_null_25, rho_null_50, rho_null_75, rho_null_95, rho_null_975
        )
   ')
    
    indb_filenames <- list.files(file.path(script_dir, 'jobs'), recursive = TRUE, pattern = '*.sqlite$')
    for(indb_filename in indb_filenames) {
        load_input_db(db, file.path(script_dir, 'jobs', indb_filename))
    }
    dbDisconnect(db)
    
    invisible()
}

load_input_db <- function(db, indb_filename) {
    dbGetPreparedQuery(db, 'ATTACH ? AS indb', data.frame(indb = indb_filename))
    dbGetQuery(db, 'INSERT INTO ccm_rho SELECT * FROM indb.ccm_rho')
    dbGetQuery(db, 'DETACH indb')
}

get_script_dir <- function()
{
    command_args <- commandArgs(trailingOnly = FALSE)
    file_arg_name <- '--file='
    tools:::file_path_as_absolute(
        dirname(
            sub(file_arg_name, '', command_args[grep(file_arg_name, command_args)])
        )
    )
}

main()
