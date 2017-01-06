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
    
    for(country_name in c(
        'Malaysia',
        'Colombia',
        'Peru',
        'Philippines',
        'Senegal',
        'Thailand',
        'Madagascar',
        'New Caledonia',
        'Paraguay',
        'South Africa',
        'Chile',
        'Japan',
        'Portugal',
        'Spain',
        'Romania',
        'France',
        'Slovenia',
        'Switzerland',
        'Germany',
        'Poland',
        'United Kingdom',
        'Denmark',
        'Latvia',
        'Sweden',
        'Norway',
        'Finland'
    )) {
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

process_remove_zeros <- function(df, env_var) {
    df[[env_var]][df$flu == 0.0] <- NA
    df$flu[df$flu == 0.0] <- NA
    df
}

process_log_flu <- function(df) {
    df$flu <- log(df$flu)
    df$flu[is.nan(df$flu)] <- NA
    df$flu[is.infinite(df$flu)] <- NA
    df
}

run_ccm_surrogates <- function(df, use_splines, cause, use_surr_cause, effect, use_surr_effect) {
    stopifnot(use_surr_cause || use_surr_effect)
    
    if(use_surr_cause) {
        cause_surrogates <- draw_surrogates(df, cause, use_splines)
    }
    if(use_surr_effect) {
        effect_surrogates <- draw_surrogates(df, effect, use_splines)
    }
    
    unlist(mclapply(1:N_SURROGATES, function(i) {
        if(use_surr_cause) {
            cause_ts <- cause_surrogates[, i]
        }
        else {
            cause_ts <- df[[cause]]
        }
        
        if(use_surr_effect) {
            effect_ts <- effect_surrogates[, i]
        }
        else {
            effect_ts <- df[[effect]]
        }
        
        run_ccm(cause_ts, effect_ts)
    }, mc.cores = N_CORES))
}

draw_surrogates <- function(df, var, use_splines) {
    if(use_splines) {
        return(draw_surrogates_splines(df$date, df[[var]]))
    }
    else {
        return(draw_surrogates_nosplines(df$week, df[[var]]))
    }
}

draw_surrogates_splines <- function(t, x) {
    x[x == 0.0] <- NA
    
    anom <- yearday_anom(t, x)
    do.call(
        cbind,
        lapply(1:N_SURROGATES, function(i) {
            I_na <- is.na(anom$anomaly)
            out <- anom$mean
            out[I_na] <- NA
            out[!I_na] <- out[!I_na] + sample(
                anom$anomaly[!I_na], sum(!I_na), replace = FALSE
            )
            return(out)
        })
    )
}

draw_surrogates_nosplines <- function(week, x) {
    anom <- weekly_anom(week, x)
    do.call(
        cbind,
        lapply(1:N_SURROGATES, function(i) {
            I_na <- is.na(anom$anomaly)
            out <- anom$mean
            out[I_na] <- NA
            out[!I_na] <- out[!I_na] + sample(
                anom$anomaly[!I_na], sum(!I_na), replace = FALSE
            )
            return(out)
        })
    )
}

run_ccm <- function(cause_ts, effect_ts) {
    block <- data.frame(effect = effect_ts, cause = cause_ts)
    lib_ccm <- c(1, length(cause_ts))
    
    out_temp <- do.call(
        rbind,
        lapply(1:8, function(E_i) {
            pred_ccm <- make_pred_nozero(effect_ts, E_i)
            ccm(
                block = block,
                E = E_i,
                lib = lib_ccm,
                pred = pred_ccm,
                lib_sizes = length(cause_ts),
                exclusion_radius = 0,
                random_libs = FALSE,
                num_sample = 1,
                tp = -1,
                lib_column = 1,
                target_column = 2
            )
        })
    )
    E_star = out_temp$E[which.max(out_temp$rho)]
    
    pred_ccm <- make_pred_nozero(effect_ts, E_star)
    ccm(
        block = block,
        E = E_star,
        lib = lib_ccm,
        pred = pred_ccm,
        lib_sizes = length(cause_ts),
        exclusion_radius = 0,
        random_libs = FALSE,
        num_sample = 1,
        tp = 0,
        lib_column = 1,
        target_column = 2
    )$rho
}

make_pred_nozero <- function(time_series, E) {
    I_zero_strings <- which(time_series == 0)
    I_zero_strings <- Reduce(
        intersect,
        lapply((0:E), function(offset) {
            I_zero_strings - offset
        })
    )
    I_zero_strings <- c(0, I_zero_strings, length(time_series))
    N_zero_strings <- length(I_zero_strings)
    lib_nozeroes <- cbind(
        I_zero_strings[1:(N_zero_strings - 1)] + 1,
        I_zero_strings[2:(N_zero_strings)]
    )
    
    lib_out <- lib_nozeroes[which(lib_nozeroes[,2] > lib_nozeroes[,1]),]
    return(lib_out)
}

yearday_anom <- function(t, x) {
    #     print('yearday_anom()')
    # t: date formatted with POSIXt
    # x: time-series values to compute seasonal mean and anomaly
    
    doy <- as.numeric(strftime(t, format = '%j'))
    I_use <- which(!is.na(x))
    
    # create time indices to use for smoothing, replicating data to "wrap around"
    doy_sm <- rep(doy[I_use], 3) + rep(c(-366, 0, 366), each = length(I_use))
    x_sm <- rep(x[I_use], 3)
    
    xsp <- smooth.spline(
        doy_sm, y = x_sm, w = NULL, spar = 0.8, cv = NA,
        all.knots = TRUE, keep.data = TRUE, df.offset = 0
    )
    
    xbar <- data.frame(t = t, doy = doy) %>%
        left_join(data.frame(doy = xsp$x, xbar = xsp$y), by='doy') %>%
        select(xbar)
    
    out <- data.frame(t = t, mean = xbar, anomaly = (x - xbar))
    names(out) <- c('date', 'mean', 'anomaly')
    
    return(out)
}

make_data_periodic <- function(df) {
    doy <- as.numeric(strftime(df$date, format = '%j'))
    df$week <- as.integer(((doy - 1) / 7) + 1)
    
    df %>% filter(week <= 52)
}

weekly_anom <- function(week, x) {
    xbar <- rep_len(sapply(1:52, function(i) {
        mean(x[week == i], na.rm = T)
    }), length.out = length(x))
    
    out <- data.frame(mean = xbar, anomaly = (x - xbar))
    names(out) <- c('mean', 'anomaly')
    
    out
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
