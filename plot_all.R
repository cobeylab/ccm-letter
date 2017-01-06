#!/usr/bin/env Rscript

library(RSQLite)
library(ggplot2)

main <- function() {
    script_dir <- get_script_dir()
    db <- dbConnect(SQLite(), file.path(script_dir, 'output_db_all.sqlite'))
    
    for(use_splines in c(FALSE, TRUE)) {
        for(remove_zeros in c(FALSE, TRUE)) {
            for(use_surr_vec in list(c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))) {
                use_surr_flu <- use_surr_vec[1]
                use_surr_env <- use_surr_vec[2]
                for(use_log_flu in c(FALSE, TRUE)) {
                    plot_model(
                        script_dir, db,
                        use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu
                    )
                }
            }
        }
    }
}

plot_model <- function(script_dir, db, use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu) {
    plot_dir <- file.path(
        script_dir, 'plots',
        sprintf(
            'us=%d-rz=%d-usf=%d-use=%d-ulf=%d',
            use_splines, remove_zeros, use_surr_flu, use_surr_env, use_log_flu
        )
    )
    if(!dir.exists(plot_dir)) {
        dir.create(plot_dir, recursive = TRUE)
    }
    
    df <- dbGetPreparedQuery(db, 
        'SELECT * FROM ccm_rho WHERE use_splines = ? AND remove_zeros = ? AND use_surr_flu = ? AND use_surr_env = ? AND use_log_flu = ? AND cause = "flu"',
        data.frame(
            use_splines = use_splines, remove_zeros = remove_zeros, use_surr_flu = use_surr_flu, use_surr_env = use_surr_env, use_log_flu = use_log_flu
        )
    )
    countries <- read.table(file.path(script_dir, 'countries.txt'), colClasses = 'character', sep = '\t')[,1]
    df$country_factor <- factor(df$country, levels = countries)
    df$effect_factor <- factor(df$effect, levels = c('AH', 'T', 'RH', 'PRCP'))
    
    p <- ggplot(data = df) +
        geom_point(aes(x = rho, y = 0, color = factor(rho_sig_95), size = 1 + rho_sig_95)) +
        scale_radius(range = c(0.5, 1)) +
        scale_color_manual(values = c('black', 'red')) +
        geom_vline(aes(xintercept = rho_null_95), size = 0.25) +
        labs(x = 'cross-map correlation', y = '') +
        facet_grid(
            country_factor ~ effect, switch = 'y',
            labeller = labeller(effect = function(effect) {
                c(AH = 'abs. humidity', RH = 'rel. humidity', T = 'temp.', PRCP = 'precip.')[effect]
            })
        ) +
        expand_limits(x = c(-0.25, 1)) +
        theme_minimal(base_size = 9) +
        theme(
            strip.text.y = element_text(angle = 180, hjust = 1),
            axis.text.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = 'none',
            panel.border = element_rect(linetype = 1, fill = NA),
            panel.spacing.x = unit(0.25, 'cm'),
            panel.spacing.y = unit(0.125, 'cm')
        )
    
    ggsave(file.path(plot_dir, 'rho_flu_cause.pdf'), p, width = 4, height = 4)
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
