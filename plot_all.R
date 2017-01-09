#!/usr/bin/env Rscript

library(RSQLite)
library(ggplot2)

main <- function() {
    script_dir <- get_script_dir()
    db <- dbConnect(SQLite(), file.path(script_dir, 'output_db_all.sqlite'))
    
    plot_model_summary(script_dir, db)
    
    for(use_splines in c(FALSE, TRUE)) {
        for(remove_zeros in c(FALSE, TRUE)) {
            for(use_surr_vec in list(c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))) {
                use_surr_flu <- use_surr_vec[1]
                use_surr_env <- use_surr_vec[2]
                for(use_log_flu in c(FALSE, TRUE)) {
                    for(flu_is_cause in c(FALSE, TRUE)) {
                        for(use_lagtest in c(FALSE, TRUE)) {
                            plot_model(
                                script_dir, db,
                                use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, flu_is_cause, use_lagtest
                            )
                        }
                    }
                }
            }
        }
    }
}

plot_model_summary <- function(script_dir, db) {
    if(!dir.exists('plots')) {
        dir.create('plots')
    }
    
    query_format <- '
        SELECT
            ccm_rho.use_surr_flu, ccm_rho.use_surr_env, ccm_rho.use_log_flu,
        SUM(rho_sig_95) AS n_sig_95,
        SUM(rho_sig_95 AND (best_lag <= 0)) AS n_sig_95_lag
        FROM ccm_rho, ccm_lagtest
        WHERE
            ccm_rho.%s = "flu" AND ccm_lagtest.%s = "flu"
            AND ccm_rho.use_splines = 1
            AND ccm_rho.remove_zeros = 1 AND ccm_lagtest.remove_zeros = 1
            AND ccm_rho.use_log_flu = ccm_lagtest.use_log_flu
            AND ccm_rho.%s = ccm_lagtest.%s
            AND ccm_rho.country = ccm_lagtest.country
        GROUP BY
            ccm_rho.use_surr_flu, ccm_rho.use_surr_env, ccm_rho.use_log_flu
        ORDER BY
            ccm_rho.use_log_flu, ccm_rho.use_surr_flu, ccm_rho.use_surr_env
    '
    
    df_flucause <- dbGetQuery(
        db,
        sprintf(query_format, 'cause', 'cause', 'effect', 'effect')
    )
    df_envcause <- dbGetQuery(
        db,
        sprintf(query_format, 'effect', 'effect', 'cause', 'cause')
    )
    stopifnot(nrow(df_flucause) == nrow(df_envcause))
    
    table_file <- file(file.path(script_dir, 'summary-table.tex'), 'w')
    format_bool <- function(x) {
        if(x) {
            return('X')
        }
        else {
            return('')
        }
    }
    for(i in 1:nrow(df_flucause)) {
        cat(
            sprintf(
                '%s & %s & %s & %d & %d & %d & %d \\\\\n',
                format_bool(df_flucause$use_surr_flu[i]),
                format_bool(df_flucause$use_surr_env[i]),
                format_bool(df_flucause$use_log_flu[i]),
                df_flucause$n_sig_95[i], df_envcause$n_sig_95[i],
                df_flucause$n_sig_95_lag[i], df_envcause$n_sig_95_lag[i]
            ),
            file = table_file
        )
    }
    close(table_file)
}

plot_model <- function(script_dir, db, use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, flu_is_cause, use_lagtest) {
    if(!dir.exists('plots')) {
        dir.create('plots')
    }
    
    if(flu_is_cause) {
        flu_cause_or_effect <- 'cause'
        var_cause_or_effect <- 'effect'
    }
    else {
        flu_cause_or_effect <- 'effect'
        var_cause_or_effect <- 'cause'
    }
    
    query <- sprintf('
        SELECT ccm_rho.cause, ccm_rho.effect, ccm_rho.country, rho_sig_95, rho, rho_null_95, best_lag, rho_best_lag
        FROM ccm_rho, ccm_lagtest
        WHERE ccm_rho.use_splines = ? AND ccm_rho.remove_zeros = ? AND ccm_rho.use_surr_flu = ?
            AND ccm_rho.use_surr_env = ? AND ccm_rho.use_log_flu = ? AND ccm_rho.%s = "flu"
            AND ccm_lagtest.remove_zeros = ? AND ccm_lagtest.use_log_flu = ? and ccm_lagtest.%s = "flu"
            AND ccm_rho.%s = ccm_lagtest.%s AND ccm_rho.country = ccm_lagtest.country
    ', flu_cause_or_effect, flu_cause_or_effect, var_cause_or_effect, var_cause_or_effect)
    
    df <- dbGetPreparedQuery(db, 
        query,
        data.frame(
            ccm_rho.use_splines = use_splines, ccm_rho.remove_zeros = remove_zeros, ccm_rho.use_surr_flu = use_surr_flu,
            ccm_rho.use_surr_env = use_surr_env, ccm_rho.use_log_flu = use_log_flu,
            ccm_lagtest.remove_zeros = remove_zeros, ccm_lagtest.use_log_flu = use_log_flu
        )
    )
    countries <- read.table(file.path(script_dir, 'countries.txt'), colClasses = 'character', sep = '\t')[,1]
    df$country_factor <- factor(df$country, levels = countries)
    if(flu_is_cause) {
        df$var_factor <- factor(df$effect, levels = c('AH', 'T', 'RH', 'PRCP'))
    }
    else {
        df$var_factor <- factor(df$cause, levels = c('AH', 'T', 'RH', 'PRCP'))
    }

    if(use_lagtest) {
        df$significant <- df$rho_sig_95 & (df$best_lag <= 0)
    }
    else {
        df$significant <- df$rho_sig_95
    }

    p <- ggplot(data = df) +
        geom_point(aes(x = rho, y = 0, color = factor(significant), size = 1 + significant)) +
        scale_radius(range = c(0.5, 1)) +
        scale_color_manual(values = c('black', 'red')) +
        geom_vline(aes(xintercept = rho_null_95), size = 0.25) +
        labs(x = 'cross-map correlation (flu drives env. variable)', y = NULL) +
        facet_grid(
            country_factor ~ var_factor, switch = 'y',
            labeller = labeller(var_factor = function(var) {
                df_sig_sum <- sapply(var, function(v) { sum(df$significant[df$var_factor == v]) })
                label <- c(
                    AH = 'abs. hum.', RH = 'rel. hum.', T = 'temp.', PRCP = 'precip.'
                )[as.character(var)]
                sprintf('%s (%d/26)', label, df_sig_sum)
            })
        ) +
        expand_limits(x = c(-0.25, 1)) +
        scale_x_continuous(breaks=c(0, 0.5, 1.0)) +
        theme_minimal(base_size = 9) +
        theme(
            strip.text.y = element_text(angle = 180, hjust = 1),
            axis.text.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_line(size = 0.25, linetype = 2, color = 'darkgray'),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = 'none',
            panel.border = element_rect(linetype = 1, fill = NA),
            panel.spacing.x = unit(0.25, 'cm'),
            panel.spacing.y = unit(0.125, 'cm'),
            axis.title.x = element_text(size = 8),
            axis.ticks.x = element_line(size = 0.25)
        )

    plot_filename <- file.path(
        script_dir, 'plots',
        sprintf(
            'us=%d-rz=%d-usf=%d-use=%d-ulf=%d-fic=%d-lt=%d.pdf',
            use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu, flu_is_cause, use_lagtest
        )
    )
    ggsave(plot_filename, p, width = 4.5, height = 4)
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
