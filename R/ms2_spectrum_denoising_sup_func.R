


get_transform_function = function(config) {
    if (config$method %in% c('fixed')) {
        ms2_spectrum_transform_fun = function(z) z
    }

    # built-in transform functions
    if (is.character(config$transform)) {
        ms2_spectrum_transform_fun = switch(config$transform, 
        log2 = function(z) log2(z+1),
        sqrt = function(z) sqrt(z),
        asinh = function(z) asinh(z),
        none = function(z) z
        )
    } else if (is.function(config$transform)) {
        ms2_spectrum_transform_fun = config$transform
    }

    return(ms2_spectrum_transform_fun)
}


prepare_ms_spectra_data_4_ms2_denoising = function(ms_data, config) {
    # check if ms_data is empty
    if (length(ms_data) == 0) {
        stop("No MS spectra found in `ms_data`.")
    }

    ms2_spectra_data = Spectra::filterMsLevel(ms_data, msLevel = 2)
    ms2_peaks_data = Spectra::peaksData(ms2_spectra_data)

    get_meta <- function(name, ms2_spectra_data, missing_value) {
        n <- length(ms2_spectra_data)
        if (!name %in% Spectra::spectraVariables(ms2_spectra_data)) {
            return(rep(missing_value, n))
        }
        value <- ms2_spectra_data[[name]]
        if (length(value) != n) stop("Metadata length mismatch: ", name)
        
        return(value)
    }

    ms2_metadata = data.frame(
        ms2_spectrum_id = ms2_spectra_data[['spectrumId']],
        ms2_acquisition_num = get_meta('acquisitionNum', ms2_spectra_data, NA),
        ms2_scan_index = get_meta('scanIndex', ms2_spectra_data, NA),
        ms2_rt = get_meta('rtime', ms2_spectra_data, NA),
        ms2_precursor_mz = get_meta('precursorMz', ms2_spectra_data, NA),
        ms2_precursor_intensity = get_meta('precursorIntensity', ms2_spectra_data, NA),
        ms2_precursor_charge = get_meta('precursorCharge', ms2_spectra_data, NA),
        ms2_peaks_count = get_meta('peaksCount', ms2_spectra_data, NA),
        ms2_tot_ion_current = get_meta('totIonCurrent', ms2_spectra_data, NA),
        ms2_base_peak_mz = get_meta('basePeakMZ', ms2_spectra_data, NA),
        ms2_base_peak_intensity = get_meta('basePeakIntensity', ms2_spectra_data, NA),
        ms2_collision_energy = get_meta('collisionEnergy', ms2_spectra_data, NA),
        ms2_ionisation_energy = get_meta('ionisationEnergy', ms2_spectra_data, NA),
        ms2_isolation_window_target_mz = get_meta('isolationWindowTargetMz', ms2_spectra_data, NA),
        ms2_ion_mobility_drift_time = get_meta('ionMobilityDriftTime', ms2_spectra_data, NA)
    )

    prepared_ms_data = structure(
        list(
            ms_data = ms_data, 
            ms2_spectra_data = ms2_spectra_data,
            ms2_peaks_data = ms2_peaks_data,
            ms2_metadata = ms2_metadata
        ), 
        class = "prepared_ms_data"
    )

    return(prepared_ms_data)

}




get_spline_denoising_info = function(current_ms2_spectra_peaks, ms2_signal_transform_func, config) {
    # preprocessing
    current_ms2_spectra_signal_intensity_distinct = dplyr::distinct(current_ms2_spectra_peaks, intensity, .keep_all = TRUE)

    # non-linear transformation of the signal intensities
    ms2_spectra_transformed_signal = sort(ms2_signal_transform_func(current_ms2_spectra_signal_intensity_distinct$intensity))

    # determine the optimal spar for the spline smoothing of the transformed signal
    spar_seq = seq(from = config$spar_start, to = config$spar_end, by = config$spar_step)

    spar_eval_data = data.frame(spar = numeric(), RMSE = numeric(), CV = numeric(), D1 = numeric(), D2 = numeric())

    # RMSE, CV, D1, D2 use 2 times denser grid for better accuracy
    x_grid = seq(from = 1, to = length(ms2_spectra_transformed_signal), length.out = length(ms2_spectra_transformed_signal)*2)

    for (spar in spar_seq) {
        # perform spline smoothing with the current spar
        spline_fit = stats::smooth.spline(
            x = seq_along(ms2_spectra_transformed_signal), 
            y = ms2_spectra_transformed_signal, 
            spar = spar, 
            cv = config$spline_cv
        )

        # calculate evaluation metrics (RMSE, CV, D1, D2)
        # RMSE
        residuals = ms2_spectra_transformed_signal - stats::predict(spline_fit, seq_along(ms2_spectra_transformed_signal))$y
        RMSE = sqrt(mean(residuals^2))

        # loocv
        if (!is.null(spline_fit$cv.crit)) {
            LOOCV = spline_fit$cv.crit
        } else {
            LOOCV = NA
        }

        # D1
        d1 = stats::predict(
            spline_fit, 
            x_grid, 
            deriv = 1)$y

        D1 = sum(abs(diff(d1)))

        # D2
        d2 = stats::predict(
            spline_fit, 
            x_grid, 
            deriv = 2)$y

        dx = diff(x_grid)
        dx = c(dx, mean(dx))

        D2 = sum((d2^2) * dx)

        spar_eval_val = data.frame(spar = spar, RMSE = RMSE, CV = LOOCV, D1 = D1, D2 = D2)

        # store the evaluation results
        spar_eval_data = dplyr::bind_rows(spar_eval_data, spar_eval_val)
    }
    
    # min max normalization
    # min_max_norm = function(x) {
    #     (x - min(x)) / (max(x) - min(x) +1e-12)
    # }
    min_max_norm <- function(x) {
        ok <- is.finite(x)
        result <- rep(NA_real_, length(x))

        if (!any(ok)) {
            return(result)
        }

        lower <- min(x[ok])
        upper <- max(x[ok])
        result[ok] <- (x[ok] - lower) / (upper - lower + 1e-12)

        return(result)
    }

    spar_eval_data_norm = data.frame(
        spar = spar_eval_data$spar,
        RMSE_norm = min_max_norm(spar_eval_data$RMSE),
        CV_norm = min_max_norm(spar_eval_data$CV),
        D1_norm = min_max_norm(spar_eval_data$D1),
        D2_norm = min_max_norm(spar_eval_data$D2)
    )

    spar_eval_data_norm = dplyr::mutate(spar_eval_data_norm, 
        eval_score = dplyr::if_else(is.na(CV_norm), 
        (config$RMSE_weight + config$CV_weight) * RMSE_norm + config$D1_weight * D1_norm + config$D2_weight * D2_norm, 
        config$RMSE_weight * RMSE_norm + config$CV_weight * CV_norm + config$D1_weight * D1_norm + config$D2_weight * D2_norm)
        )

    valid_idx <- which(is.finite(spar_eval_data_norm$eval_score))

    if (length(valid_idx) == 0L) {
        stop("No finite spline candidate score.", call. = FALSE)
    }

    # determine the optimal spar
    optimal_idx <- valid_idx[which.min(spar_eval_data_norm$eval_score[valid_idx])]
    optimal_spar <- spar_eval_data_norm$spar[optimal_idx]

    # fit the optimal smoothing spline
    spline_fit_optimal = stats::smooth.spline(
        x = seq_along(ms2_spectra_transformed_signal), 
        y = ms2_spectra_transformed_signal,  
        spar = optimal_spar, 
        cv = config$spline_cv
    )

    # signal index
    d1 = stats::predict(
        spline_fit_optimal, 
        x_grid, 
        deriv = 1)$y

    # remove the top peaks based on the configuration
    top_peaks_idx = length(ms2_spectra_transformed_signal) - config$top_peaks_to_remove

    keep_idx = which(x_grid <= top_peaks_idx & is.finite(d1))
    x_grid_search = x_grid[keep_idx]
    d1_search = d1[keep_idx]

    # get the noise threshold
    threshold_rank = x_grid_search[which.min(abs(d1_search))]
    threshold_idx = as.integer(floor(threshold_rank))
    noise_threshold = sort(current_ms2_spectra_signal_intensity_distinct$intensity)[threshold_idx]

    # output
    return(
        list(
            spar_eval_data = spar_eval_data,
            spar_eval_data_norm = spar_eval_data_norm,
            optimal_spar = optimal_spar,
            spline_fit_optimal = spline_fit_optimal, 
            spline_x = seq_along(ms2_spectra_transformed_signal), 
            spline_y = ms2_spectra_transformed_signal, 
            threshold_idx = threshold_idx,
            noise_threshold = noise_threshold
            )
        )
}




get_segmented_denoising_info = function(current_ms2_spectra_peaks, ms2_signal_transform_func, config) {
    # preprocessing
    current_ms2_spectra_signal_intensity_distinct = dplyr::distinct(current_ms2_spectra_peaks, intensity, .keep_all = TRUE)

    # non-linear transformation of the signal intensities
    ms2_spectra_transformed_signal = sort(ms2_signal_transform_func(current_ms2_spectra_signal_intensity_distinct$intensity))

    x_idx = seq_along(ms2_spectra_transformed_signal)

    # lm fit 
    lm_fit = stats::lm(ms2_spectra_transformed_signal ~ x_idx)

    # segmented regression to find break point
    lm_fit_seg = segmented::segmented(lm_fit, seg.Z = ~x_idx, psi = list(x_idx = length(ms2_spectra_transformed_signal)/3))

    # get the noise threshold
    break_point = lm_fit_seg$psi[,"Est."]
    threshold_idx = round(break_point)
    noise_threshold = sort(current_ms2_spectra_signal_intensity_distinct$intensity)[threshold_idx]

    # linear and non_linear fitting
    # linear fitting
    y_linear = ms2_spectra_transformed_signal[1:threshold_idx]
    x_linear = x_idx[1:threshold_idx]

    lm_fit_linear = stats::lm(y_linear ~ x_linear)

    # non-linear fitting
    segmented_nonlinear_fit_func = config$segmented_nonlinear_fit_function

    y_non_linear = ms2_spectra_transformed_signal[(threshold_idx+1):length(ms2_spectra_transformed_signal)]
    x_non_linear = segmented_nonlinear_fit_func(
        x_idx[(threshold_idx+1):length(ms2_spectra_transformed_signal)]
    )
    lm_fit_non_linear = stats::lm(y_non_linear ~ x_non_linear)

    # output
    return(
        list(
            lm_fit = lm_fit, 
            lm_fit_linear = lm_fit_linear,
            lm_fit_non_linear = lm_fit_non_linear, 
            segmented_x = x_idx, 
            segmented_y = ms2_spectra_transformed_signal, 
            threshold_idx = threshold_idx,
            noise_threshold = noise_threshold
            )
        )
}



get_quantile_denoising_info = function(current_ms2_spectra_peaks, ms2_signal_transform_func, config) {
    # preprocessing
    current_ms2_spectra_signal_intensity_distinct = dplyr::distinct(current_ms2_spectra_peaks, intensity, .keep_all = TRUE)

    # non-linear transformation of the signal intensities
    ms2_spectra_transformed_signal = sort(ms2_signal_transform_func(current_ms2_spectra_signal_intensity_distinct$intensity))
    
    # noise threshold
    quantile_threshold = stats::quantile(ms2_spectra_transformed_signal, config$quantile_prob)

    noise_index = max(which(ms2_spectra_transformed_signal <= quantile_threshold))

    noise_threshold = sort(current_ms2_spectra_signal_intensity_distinct$intensity)[noise_index]

    # output
    return(
        list(
            noise_threshold = noise_threshold
        )
    )
}

get_fixed_denoising_info = function(current_ms2_spectra_peaks, ms2_signal_transform_func, config) {
    # noise threshold
    noise_threshold = config$fixed_value

    # output
    return(
        list(
            noise_threshold = noise_threshold
        )
    )
}






calculate_ms2_noise_threshold = function(current_ms2_spectra_peaks, ms2_signal_transform_func, config) {
    # check the number of top peaks
    if ((config$method == "adaptive") && 
        nrow(current_ms2_spectra_peaks) < config$top_peaks_for_pattern + config$top_peaks_to_remove + 4L) {
            return(list(
                denoising_method = NA,
                noise_threshold = 0, 
                denoising_status = 'fail'
            )
        )
    } else if (config$method == "spline" && 
        nrow(current_ms2_spectra_peaks) < config$top_peaks_to_remove + 4L) {
            return(list(
                denoising_method = NA,
                noise_threshold = 0, 
                denoising_status = 'fail'
            )
        )
    } else if (config$method == "segmented" && nrow(current_ms2_spectra_peaks) < 10) {
            return(list(
                denoising_method = NA,
                noise_threshold = 0, 
                denoising_status = 'fail'
            )
        )
    } else {
        # determine the denoising method
        if(config$method == "adaptive") {
            # get the top n signal intensities
            top_n_signal = dplyr::arrange(current_ms2_spectra_peaks, dplyr::desc(intensity))
            top_n_signal = head(top_n_signal, config$top_peaks_for_pattern)

            top_n_signal_mz = sort(top_n_signal$mz)

            # check if the top n signals are isotope spacing 
            is_evenly_spaced <- function(x, tol) {
                x <- sort(x)
                gaps = diff(x)
                
                if (any(gaps <= 0)) {
                    return(FALSE)
                }

                reference_gap = stats::median(gaps)

                if (all(abs(gaps - reference_gap) / reference_gap <= tol)) {
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            }

            if (is_evenly_spaced(top_n_signal_mz, config$top_peaks_evenly_spaced_tol)) {
                # use spline denoising method 
                denoising_method = "spline"
                spline_denoising_info = get_spline_denoising_info(
                    current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
                    ms2_signal_transform_func = ms2_signal_transform_func, 
                    config = config
                    )
                denoising_status = 'success'
                noise_threshold = spline_denoising_info$noise_threshold
            } else {
                # use segmented denoising method
                denoising_method = "segmented"
                segmented_denoising_info = get_segmented_denoising_info(
                    current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
                    ms2_signal_transform_func = ms2_signal_transform_func, 
                    config = config
                    )
                denoising_status = 'success'
                noise_threshold = segmented_denoising_info$noise_threshold
                }        

        } else if (config$method == "spline") {
            denoising_method = "spline"
            spline_denoising_info = get_spline_denoising_info(
                current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
                ms2_signal_transform_func = ms2_signal_transform_func, 
                config = config
            )
            denoising_status = 'success'
            noise_threshold = spline_denoising_info$noise_threshold

        } else if (config$method == "segmented") {
            denoising_method = "segmented"
            segmented_denoising_info = get_segmented_denoising_info(
                current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
                ms2_signal_transform_func = ms2_signal_transform_func, 
                config = config
            )
            denoising_status = 'success'
            noise_threshold = segmented_denoising_info$noise_threshold

        } else if (config$method == "quantile") {
            denoising_method = "quantile"
            quantile_denoising_info = get_quantile_denoising_info(
                current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
                ms2_signal_transform_func = ms2_signal_transform_func, 
                config = config
            )
            denoising_status = 'success'
            noise_threshold = quantile_denoising_info$noise_threshold

        } else if (config$method == "fixed") {
            denoising_method = "fixed"
            fixed_denoising_info = get_fixed_denoising_info(
                current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
                ms2_signal_transform_func = ms2_signal_transform_func, 
                config = config
            )
            denoising_status = 'success'
            noise_threshold = fixed_denoising_info$noise_threshold

        } else {
            stop("Invalid denoising method specified in the config.")
        }

        return(
            list(
                denoising_method = denoising_method,
                noise_threshold = noise_threshold, 
                denoising_status = denoising_status
                )
            )
    }
}







do_ms2_spectra_denoising = function(ms_data_prepared, ms2_signal_transform_func, config) {
    # check class
    if (!inherits(ms_data_prepared, "prepared_ms_data")) {
        stop("`ms_data_prepared` must be of class 'prepared_ms_data'.")
    }

    # initialize result
    ms2_length = length(ms_data_prepared$ms2_spectra_data)

    denoising_result = data.frame(
        ms2_spectrum_id = character(ms2_length),
        ms2_noise_threshold = rep(NA_real_, ms2_length),
        denoising_method = character(ms2_length),
        denoising_status = character(ms2_length),
        stringsAsFactors = FALSE
    )
    
    # check if the number of MS2 spectra matches the number of MS2 peaks
    if (ms2_length != length(ms_data_prepared$ms2_peaks_data)) {
        stop("Mismatch between the number of MS2 spectra and MS2 peaks.")
    }

    # progress bar
    pb <- utils::txtProgressBar(min = 0, max = ms2_length, style = 3)
    on.exit(close(pb))

    # start denoising loop
    for (i in seq_len(ms2_length)) {
        # update progress bar
        utils::setTxtProgressBar(pb, i)

        # extract the current MS2 spectrum and its corresponding peaks
        current_ms2_spectra = ms_data_prepared$ms2_spectra_data[i]
        current_ms2_spectra_peaks = as.data.frame(ms_data_prepared$ms2_peaks_data[[i]])

        # validate ids
        if (current_ms2_spectra$spectrumId != ms_data_prepared$ms2_metadata$ms2_spectrum_id[i]) {
            stop(paste(
                "Mismatch between the current MS2 spectrumId and the metadata.", 
                "Expected: ", ms_data_prepared$ms2_metadata$ms2_spectrum_id[i], 
                " Got: ", current_ms2_spectra$spectrumId
                )
            )
        }
        if (current_ms2_spectra$acquisitionNum != ms_data_prepared$ms2_metadata$ms2_acquisition_num[i]) {
            stop(paste(
                "Mismatch between the current MS2 spectrum acquisitionNum and the metadata.", 
                "Expected: ", ms_data_prepared$ms2_metadata$ms2_acquisition_num[i], 
                " Got: ", current_ms2_spectra$acquisitionNum
                )
            )
        }
        if (current_ms2_spectra$scanIndex != ms_data_prepared$ms2_metadata$ms2_scan_index[i]) {
            stop(paste(
                "Mismatch between the current MS2 spectrum scanIndex and the metadata.", 
                 "Expected: ", ms_data_prepared$ms2_metadata$ms2_scan_index[i], 
                 " Got: ", current_ms2_spectra$scanIndex
                )
            )
        }

        # required columns in the peaks data
        required_columns = c("mz", "intensity")

        if (!all(required_columns %in% colnames(current_ms2_spectra_peaks))) {
            stop(paste(
                "The MS2 peaks data must contain the following columns: 'mz' and 'intensity'.", 
                "Current spectrumId: ", current_ms2_spectra$spectrumId
                )
            )
        }

        denoising_info = calculate_ms2_noise_threshold(
            current_ms2_spectra_peaks = current_ms2_spectra_peaks, 
            ms2_signal_transform_func = ms2_signal_transform_func, 
            config = config)

        noise_threshold = denoising_info$noise_threshold

        denoising_result$ms2_spectrum_id[i] = as.character(current_ms2_spectra$spectrumId)
        denoising_result$ms2_noise_threshold[i] = as.numeric(noise_threshold)
        denoising_result$denoising_method[i] = denoising_info$denoising_method
        denoising_result$denoising_status[i] = denoising_info$denoising_status
    }

    if (anyDuplicated(denoising_result$ms2_spectrum_id)) {
        stop("Duplicate spectrum IDs in `denoising_result`.")
    }

    # perform denoising on the ms data and output the denoised ms data
    ms_data_noise_filter = function(x, spectrumMsLevel, spectrumId, threshold_lookup, ...) {
        if (!isTRUE(spectrumMsLevel == 2L) || is.null(dim(x)) || nrow(x) == 0L || !"intensity" %in% colnames(x)) {
            return(x)
        }

        if (length(spectrumId) != 1L || is.na(spectrumId)) {
            return(x)
        }

        idx_row = dplyr::filter(threshold_lookup, 
            as.character(.data$ms2_spectrum_id) == as.character(.env$spectrumId))
        
        if (nrow(idx_row) == 0L) {
            return(x) 
        }

        if (nrow(idx_row) > 1L) {
            stop(paste(
                "Multiple noise thresholds found for spectrum:",
                spectrumId
            ))
        }

        current_thres = idx_row$ms2_noise_threshold

        if (!is.finite(current_thres)) {
            return(x)
        }

        intensity = x[, "intensity"]
        keep = is.finite(intensity) & intensity > current_thres

        return(x[keep, , drop = FALSE])
    }

    sps = Spectra::addProcessing(
        object = ms_data_prepared$ms_data,
        FUN = ms_data_noise_filter,
        threshold_lookup = denoising_result,
        spectraVariables = c("msLevel", "spectrumId")
    )

    # recalculate the ms2 TIC and peakscount for the denoised ms2 spectra
    sps_ms2 = Spectra::filterMsLevel(sps, msLevel = 2L)
    sps_ms2_length = length(sps_ms2)

    sps_ms2_idx <- which(Spectra::msLevel(sps) == 2L)

    if (sps_ms2_length != length(Spectra::filterMsLevel(ms_data_prepared$ms_data, msLevel = 2L))) {
        stop("Mismatch in the number of MS2 spectra after denoising.")
    }
    if (sps_ms2_length != nrow(denoising_result)) {
        stop("Mismatch in the number of MS2 spectra and denoising_results.")
    }
    if (!identical(as.character(sps$spectrumId[sps_ms2_idx]), as.character(denoising_result$ms2_spectrum_id))) {
        stop("Mismatch in MS2 spectrum IDs or their order.")
    }

    success_denoised_ms2_id = dplyr::filter(denoising_result, denoising_status == "success")$ms2_spectrum_id
    success_denoised_ms2_idx = which(as.character(sps[["spectrumId"]]) %in% as.character(success_denoised_ms2_id))
    success_denoised_ms2_peaks_data <- Spectra::peaksData(sps[success_denoised_ms2_idx])

    # recalculate peaksCount
    denoised_peaks_count <- vapply(
        success_denoised_ms2_peaks_data,
        function(p) nrow(p),
        integer(1)
        )
    # recalculate TIC; the sum of intensities for empty peak matrices is 0
    denoised_tic <- vapply(
        success_denoised_ms2_peaks_data,
        function(p) sum(p[, "intensity"]),
        numeric(1)
        )

    # update peaksCount and TIC
    updated_peaks_count <- sps[["peaksCount"]]
    updated_tic <- sps[["totIonCurrent"]]

    updated_peaks_count[success_denoised_ms2_idx] <- denoised_peaks_count
    updated_tic[success_denoised_ms2_idx] <- denoised_tic

    sps[["peaksCount"]] <- updated_peaks_count
    sps[["totIonCurrent"]] <- updated_tic


    return(
        list(denoising_result = denoising_result, 
            denoised_ms_data = sps
            )
        )
}




plot_ms2_spectra = function(ms2_spectrum_id, ms_data, denoised_data) {
    ms2_data = Spectra::filterMsLevel(ms_data, msLevel = 2)

    ms2_ids = ms2_data[["spectrumId"]]
    ms2_idx = which(ms2_ids == ms2_spectrum_id)
    spectra_data_raw = as.data.frame(Spectra::peaksData(ms2_data)[[ms2_idx]])

    denoising_detail = dplyr::filter(denoised_data$denoising_result, .data$ms2_spectrum_id == .env$ms2_spectrum_id)

    # ms2 spectrum plot
    if (denoising_detail$denoising_status == "success") {
        # denoised spectra plot
        spectra_plot = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) + 
            ggplot2::geom_segment(
                ggplot2::aes(xend = mz, yend = 0, colour = ifelse(intensity <= denoising_detail$ms2_noise_threshold, "#E0BA48", "black")), 
                linewidth = 0.2) + 
            ggplot2::scale_colour_identity() +
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.05))) +
            ggplot2::labs(
                title = paste("MS2 spectrum: ", ms2_spectrum_id, "(Denoising status: ", denoising_detail$denoising_status, ")"), 
                x = "m/z", y = "Intensity") +
            ggplot2::theme_classic()
    } else {
        spectra_plot = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) + 
            ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0), colour = "black", linewidth = 0.2) + 
            ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.05))) +
            ggplot2::labs(
                title = paste("MS2 spectrum: ", ms2_spectrum_id, "(Denoising status: ", denoising_detail$denoising_status, ")"), 
                x = "m/z", y = "Intensity") +
            ggplot2::theme_classic()
    }
    return(spectra_plot)
}




plot_ms2_denoising_info = function(ms2_spectrum_id, ms_data, denoised_data, ms2_signal_transform_func) {
    ms2_data = Spectra::filterMsLevel(ms_data, msLevel = 2)

    ms2_ids = ms2_data[["spectrumId"]]
    ms2_idx = which(ms2_ids == ms2_spectrum_id)
    spectra_data_raw = as.data.frame(Spectra::peaksData(ms2_data)[[ms2_idx]])

    denoising_detail = dplyr::filter(denoised_data$denoising_result, .data$ms2_spectrum_id == .env$ms2_spectrum_id)

    if (denoising_detail$denoising_method == "spline") {
        spline_denoising_info = get_spline_denoising_info(
            current_ms2_spectra_peaks = spectra_data_raw, 
            ms2_signal_transform_func = ms2_signal_transform_func, 
            config = denoised_data$config)

        # spar eval data plot
        spar_eval_data = spline_denoising_info$spar_eval_data
        optimal_spar = spline_denoising_info$optimal_spar

        spar_eval_data = dplyr::arrange(spar_eval_data, spar)

        metric_plots = lapply(c("RMSE", "CV", "D1", "D2"), function(metric) {
            metric_plot_data = dplyr::select(spar_eval_data, spar, !!rlang::sym(metric))
            metric_plot_data = dplyr::filter(metric_plot_data, is.finite(!!rlang::sym(metric)))
            
            ggplot2::ggplot(metric_plot_data, ggplot2::aes(x = spar, y = !!rlang::sym(metric))) + 
            ggplot2::geom_point(colour = "black", size = 0.7) + 
            ggplot2::geom_vline(xintercept = optimal_spar, colour = "#D97732", linetype = "dashed", linewidth = 0.4) + 
            ggplot2::scale_x_continuous(
                limits = range(spar_eval_data$spar),
                expand = ggplot2::expansion(mult = c(0.04, 0.04))
                ) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.03, 0.05))
                ) +
            ggplot2::labs(
                title = paste(metric),
                x = "spar",
                y = paste(metric)
                ) +
            ggplot2::theme_classic(base_size = 12) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
                axis.title = ggplot2::element_text(size = 12),
                axis.text = ggplot2::element_text(colour = "black")
                )
            }
        )
        # combine
        spline_spar_metric_plots = 
            patchwork::wrap_plots(metric_plots, nrow = 1) + 
            patchwork::plot_annotation(title = ms2_spectrum_id)

        # optimal spline plot
        threshold_idx = spline_denoising_info$threshold_idx

        optimal_spline_y = stats::predict(spline_denoising_info$spline_fit_optimal, spline_denoising_info$spline_x)$y

        orginal_spline_x_y_data = data.frame(
            x = spline_denoising_info$spline_x,
            y = spline_denoising_info$spline_y
            )

        optimal_spline_plot_data = data.frame(
            x = spline_denoising_info$spline_x, 
            y = optimal_spline_y
            )

        optimal_spline_plot = ggplot2::ggplot() + 
            ggplot2::geom_line(data = optimal_spline_plot_data[optimal_spline_plot_data$x <= threshold_idx, ,drop = FALSE], 
                ggplot2::aes(x = x, y = y), colour = "#A8B1C1", linewidth = 2) + 
            ggplot2::geom_line(data = optimal_spline_plot_data[optimal_spline_plot_data$x >= threshold_idx, , drop = FALSE],
                ggplot2::aes(x = x, y = y), color = "#8DBFDF", linewidth = 2) + 
            ggplot2::geom_vline(
                xintercept = threshold_idx, 
                linetype = "dashed", 
                color = "#E0873F", 
                linewidth = 0.5) + 
            ggplot2::geom_point(data = orginal_spline_x_y_data, ggplot2::aes(x = x, y = y), color = "black", size = 1) + 
            ggplot2::labs(title = ms2_spectrum_id, x = "Index", y = "Transformed intensity") +
            ggplot2::theme_classic(base_size = 18) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                axis.text = ggplot2::element_text(color = "black"),
                axis.title = ggplot2::element_text(color = "black")
            )

        return(
            list(
                spar_metric_plots = spline_spar_metric_plots,
                optimal_spline_plot = optimal_spline_plot
            )
        )

    } else if (denoising_detail$denoising_method == "segmented") {
        segmented_denoising_info = get_segmented_denoising_info(
            current_ms2_spectra_peaks = spectra_data_raw, 
            ms2_signal_transform_func = ms2_signal_transform_func, 
            config = denoised_data$config)

        # points data
        orginal_segmented_x_y_data = data.frame(
            x = segmented_denoising_info$segmented_x,
            y = segmented_denoising_info$segmented_y
            )

        # get the non-linear transform func
        segmented_nonlinear_fit_func = denoised_data$config$segmented_nonlinear_fit_function

        # non-linear data
        non_linear_x = segmented_denoising_info$segmented_x[(segmented_denoising_info$threshold_idx+1):length(segmented_denoising_info$segmented_x)]
        non_linear_x_transformed = segmented_nonlinear_fit_func(
            segmented_denoising_info$segmented_x[(segmented_denoising_info$threshold_idx+1):length(segmented_denoising_info$segmented_x)]
            )
        segmented_regression_nonlinear_plot_data = data.frame(
            x = non_linear_x, 
            y = stats::predict(segmented_denoising_info$lm_fit_non_linear, newdata = data.frame(x_non_linear_transformed = non_linear_x_transformed))
        )

        # linear data
        linear_x = segmented_denoising_info$segmented_x[1:segmented_denoising_info$threshold_idx]
        segmented_regression_linear_plot_data = data.frame(
            x = linear_x,
            y = stats::predict(segmented_denoising_info$lm_fit_linear, newdata = data.frame(x_linear = linear_x))
        )

        # initial linear regression plot data
        initial_linear_regression_plot_data = data.frame(
            x = segmented_denoising_info$segmented_x,
            y = stats::predict(segmented_denoising_info$lm_fit, newdata = data.frame(x = segmented_denoising_info$segmented_x))
        )

        segmented_regression_plot = ggplot2::ggplot() +
            ggplot2::geom_line(data = initial_linear_regression_plot_data, ggplot2::aes(x = x, y = y), color = "#D5E5EA", linewidth = 2) + 
            ggplot2::geom_line(data = segmented_regression_linear_plot_data, ggplot2::aes(x = x, y = y), color = "#8C96A8", linewidth = 2) + 
            ggplot2::geom_line(data = segmented_regression_nonlinear_plot_data, ggplot2::aes(x = x, y = y), color = "#CD8989", linewidth = 2) + 
            ggplot2::geom_vline(
                        xintercept = segmented_denoising_info$threshold_idx, 
                        linetype = "dashed", 
                        color = "#E0873F", 
                        linewidth = 0.5) + 
            ggplot2::geom_point(data = orginal_segmented_x_y_data, ggplot2::aes(x = x, y = y), color = "black", size = 1) + 
            ggplot2::labs(title = ms2_spectrum_id, x = "Index", y = "Transformed intensity") +
            ggplot2::theme_classic(base_size = 18) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = 12, hjust = 0.5),
                axis.text = ggplot2::element_text(color = "black"),
                axis.title = ggplot2::element_text(color = "black")
                )

        return(
            list(
                segmented_regression_plot = segmented_regression_plot
            )
        )
    } else{
        return(NULL)
    }
}









