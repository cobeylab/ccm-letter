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
    E_star = identify_E_star(cause_ts, effect_ts, block, lib_ccm)
    
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

run_ccm_lagtest <- function(cause_ts, effect_ts, lag_range) {
    block <- data.frame(effect = effect_ts, cause = cause_ts)
    lib_ccm <- c(1, length(cause_ts))
    E_star = identify_E_star(cause_ts, effect_ts, block, lib_ccm)
    
    ccm_by_lag <- do.call(
        rbind,
        lapply(lag_range, function(lag) {
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
                tp = lag,
                lib_column = 1,
                target_column = 2
            )
        })
    )
    
    best_lag_index <- which.max(ccm_by_lag$rho)
    list(
        best_lag = lag_range[best_lag_index],
        rho_best_lag = ccm_by_lag$rho[best_lag_index]
    )
}

identify_E_star <- function(cause_ts, effect_ts, block, lib_ccm) {
    result <- do.call(
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
    result$E[which.max(result$rho)]
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

load_countries <- function() {
    read.table(file.path(get_script_dir(), 'countries.txt'), colClasses = 'character', sep = '\t')[,1]
}
