#!/usr/bin/env Rscript

library(rEDM)
library(dplyr)
library(tidyr)
library(parallel)
library(RSQLite)

LAG_RANGE <- -20:20

main <- function() {
    df_all <- read.csv('pnas.1607747113.sd01.txt') %>% filter(year >= 1996)
    script_dir <- get_script_dir()
    
    df_in <- read.csv(file.path(script_dir, 'pnas.1607747113.sd01.txt')) %>% filter(year >= 1996)
    
    db_out <- dbConnect(SQLite(), file.path(script_dir, 'output_db_all.sqlite'))
    dbGetQuery(db_out, 'DROP TABLE IF EXISTS ccm_lagtest')
    dbGetQuery(db_out, '
       CREATE TABLE ccm_lagtest (
           remove_zeros, use_log_flu, cause, effect, country,
           best_lag,
           rho_best_lag
       )
    ')
    
    for(remove_zeros in c(FALSE, TRUE)) {
        for(use_log_flu in c(FALSE, TRUE)) {
            for(env_var in c('T', 'AH', 'RH', 'PRCP')) {
                for(flu_is_cause in c(FALSE, TRUE)) {
                    for(country_name in load_countries()) {
                        df_country <- df_in %>% filter(country == country_name) %>%
                            select(-country) %>%
                            spread(variable, value) %>%
                            mutate(date = ISOdate(year, month, day)) %>%
                            select(-year, -month, -day) %>%
                            select(date, flu, everything())
                        
                        run_lag_test(
                            df_country, db_out,
                            remove_zeros, use_log_flu, env_var, flu_is_cause,
                            country_name
                        )
                    }
                }
            }
        }
    }
}

run_lag_test <- function(df_in, db_out, remove_zeros, use_log_flu, env_var, flu_is_cause, country) {
    cat(sprintf(
        'Running remove_zeros=%d, use_log_flu=%d, env_var=%s, flu_is_cause=%d, country=%s...\n',
        remove_zeros, use_log_flu, env_var, flu_is_cause, country
    ))
    
    if(remove_zeros) {
        df_in <- process_remove_zeros(df_in, env_var)
    }
    
    if(use_log_flu) {
        df_in <- process_log_flu(df_in)
    }
    
    if(flu_is_cause) {
        cause <- 'flu'
        effect <- env_var
        result <- run_ccm_lagtest(df_in$flu, df_in[[env_var]], LAG_RANGE)
    }
    else {
        cause <- env_var
        effect <- 'flu'
        result <- run_ccm_lagtest(df_in$flu, df_in[[env_var]], LAG_RANGE)
    }
    
    dbGetPreparedQuery(db_out,
       'INSERT INTO ccm_lagtest VALUES (?,?,?,?,?,?,?)',
       data.frame(
           remove_zeros = remove_zeros,
           use_log_flu = use_log_flu,
           cause = cause, effect = effect, country = country,
           best_lag = result$best_lag,
           rho_best_lag = result$rho_best_lag
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

source(file.path(get_script_dir(), 'shared.R'))
main()
