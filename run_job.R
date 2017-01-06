#!/usr/bin/env Rscript

library(rEDM)
library(dplyr)
library(tidyr)
library(RSQLite)
library(parallel)

N_CORES <- 1
N_SURROGATES <- 500

main <- function() {
    job_spec <- readRDS('job_spec.Rds')
    script_dir <- get_script_dir()
    
    df_in <- read.csv(file.path(script_dir, 'pnas.1607747113.sd01.txt')) %>% filter(year >= 1996)
    
    db_out <- dbConnect(SQLite(), 'output_db.sqlite')
    dbGetQuery(db_out, '
        CREATE TABLE ccm_rho (
            use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, cause, effect, country,
            rho_sig_95,
            rho, rho_null_025, rho_null_05, rho_null_25, rho_null_50, rho_null_75, rho_null_95, rho_null_975
        )
    ')
    
    countries <- load_countries()
    for(country_name in countries) {
        df_country <- df_in %>% filter(country == country_name) %>%
            select(-country) %>%
            spread(variable, value) %>%
            mutate(date = ISOdate(year, month, day)) %>%
            select(-year, -month, -day) %>%
            select(date, flu, everything())
        
        run_country(
            df_country, db_out,
            job_spec$use_splines, job_spec$remove_zeros, job_spec$use_surr_flu, job_spec$use_surr_env, job_spec$use_log_flu, job_spec$env_var, job_spec$flu_is_cause,
            country_name
        )
    }
}

run_country <- function(df_in, db_out, use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, env_var, flu_is_cause, country) {
    cat(sprintf('Running %s...\n', country))
    
    if(!use_splines) {
        df_in <- make_data_periodic(df_in)
    }
    
    if(remove_zeros) {
        df_in <- process_remove_zeros(df_in, env_var)
    }
    
    if(use_log_flu) {
        df_in <- process_log_flu(df_in)
    }
    
    if(flu_is_cause) {
        cause <- 'flu'
        effect <- env_var
        
        rhos_null <- run_ccm_surrogates(df_in, use_splines, 'flu', use_surr_flu, env_var, use_surr_env)
        rho <- run_ccm(df_in$flu, df_in[[env_var]])
    }
    else {
        cause <- env_var
        effect <- 'flu'
        
        rhos_null <- run_ccm_surrogates(df_in, use_splines, env_var, use_surr_env, 'flu', use_surr_flu)
        rho <- run_ccm(df_in[[env_var]], df_in$flu)
    }
    
    dbGetPreparedQuery(db_out,
        'INSERT INTO ccm_rho VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
        data.frame(
            use_splines = use_splines, remove_zeros = remove_zeros,
            use_surr_flu = use_surr_flu, use_surr_env = use_surr_env,
            use_log_flu = use_log_flu,
            cause = cause, effect = effect, country = country,
            rho_sig = rho > quantile(rhos_null, 0.95),
            rho = rho,
            rho_null_025 = quantile(rhos_null, 0.025),
            rho_null_05 = quantile(rhos_null, 0.05),
            rho_null_25 = quantile(rhos_null, 0.25),
            rho_null_50 = quantile(rhos_null, 0.50),
            rho_null_75 = quantile(rhos_null, 0.75),
            rho_null_95 = quantile(rhos_null, 0.95),
            rho_null_975 = quantile(rhos_null, 0.975)
        )
    )
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

SCRIPT_DIR <- get_script_dir()
source(file.path(get_script_dir, 'shared.R'))
main()
