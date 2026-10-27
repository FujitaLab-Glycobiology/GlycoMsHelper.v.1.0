ms_data_checker = function(ms_data) {
    ms1_data = Spectra::filterMsLevel(ms_data, 1)
    ms2_data = Spectra::filterMsLevel(ms_data, 2)

    if (length(ms1_data) == 0) {
        stop("The provided file contains no MS1 spectra.", call. = FALSE)
    }

    if (length(ms2_data) == 0) {
        stop("The provided file contains no MS2 spectra.", call. = FALSE)
    }

    charges = Spectra::precursorCharge(ms2_data)

    if (all(is.na(charges)) || all(charges == 0)) {
        stop("Error: The .mzML file lacks 'precursorCharge' information.
            Downstream analysis requires charge states.")
    } else {
        message("Check passed: MS2 charge state information is present.")
    }
}





find_jump_point_by_mean_n_sd = function(data, n) {
    # Sort the data and calculate the differences between consecutive values
    # the jump point is identified as the first value where the difference exceeds the threshold
    # threshold is defined as mean + n * standard deviation of the differences

    data_sorted = data.frame(value = sort(data))
    data_sorted = dplyr::mutate(data_sorted, data_diff = c(0, diff(value)))    

    mean_diff_val = mean(data_sorted$data_diff)
    sd_val = stats::sd(data_sorted$data_diff)

    threshold = mean_diff_val + n * sd_val

    data_possible_jump_points = dplyr::filter(data_sorted, data_diff > threshold)

    if (nrow(data_possible_jump_points) == 0) {
        stop(paste0("Jump point not found in current mean + n * sd threshold (n = ", n, ")."), call. = FALSE)
    }

    data_jump_point_value = data_possible_jump_points$value[which.max(data_possible_jump_points$value)]

    return(data_jump_point_value)
}



do_quantile_filtering = function(data, quantile_prob) {
    threshold_value = as.numeric(stats::quantile(data, probs = c(quantile_prob), na.rm = TRUE))
    return(threshold_value)
}














do_ms_data_preprocessing = function(ms_data, config) {
    # Each method returns an inclusive lower and upper bound.
    threshold_functions = list(
        mean_sd = function(data, threshold) {
            c(find_jump_point_by_mean_n_sd(data, n = threshold), Inf)
        },
        quantile = function(data, threshold) {
            do_quantile_filtering(data, quantile_prob = threshold)
        },
        start_end = function(data, threshold) {
            threshold
        }
    )

    # Subset the original object once, preserving scan order and other MS levels.
    keep = rep(TRUE, length(ms_data))
    ms_levels = ms_data[["msLevel"]]
    spectrum_ids = ms_data[["spectrumId"]]
    filter_records = list()

    for (ms_level in c(1L, 2L)) {
        level_name = paste0("ms", ms_level)
        level_config = config[[level_name]]
        level_indices = which(ms_levels == ms_level)

        level_data = ms_data[level_indices]
        level_keep = rep(TRUE, length(level_indices))

        filter_records[[level_name]] = list()

        for (variable_name in names(level_config)) {
            method_info = level_config[[variable_name]]
            method = method_info[["method"]]
            threshold_function = threshold_functions[[method]]

            values = level_data[[variable_name]]
            finite_values = is.finite(values)

            variable_keep = rep(FALSE, length(values))

            if (any(finite_values)) {
                bounds = tryCatch(
                    threshold_function(values[finite_values], method_info[["threshold"]]),
                    error = function(e) {
                        stop("Filtering ", level_name, ".", variable_name, " with '", method, "': ", conditionMessage(e),
                            call. = FALSE)
                    }
                )
                
                if (!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds) || bounds[1] > bounds[2]) {
                    stop("Invalid filtering bounds for ", level_name, ".", variable_name, ".", call. = FALSE)
                }

                variable_keep[finite_values] = values[finite_values] >= bounds[1] & values[finite_values] <= bounds[2]
            }

            filter_records[[level_name]][[variable_name]] = data.frame(
                spectrum_index = level_indices,
                spectrumId = spectrum_ids[level_indices],
                variable_name = rep(variable_name, length(values)),
                values = values,
                status = ifelse(variable_keep, "keep", "drop"),
                stringsAsFactors = FALSE
            )

            level_keep = level_keep & variable_keep
        }

        keep[level_indices] = level_keep
    }

    return(list(
        processed_ms_data = ms_data[keep], 
        filter_records = filter_records
        )
    )

}






















plot_processed_ms_data_variables = function(record_df) {
    do_plot = function(data_df, variable_name, title, x_label) {
        data_df = data_df[order(data_df$status == "keep"), , drop = FALSE]

        plot = ggplot2::ggplot(data_df, ggplot2::aes(x = spectrum_index, y = values, colour = status)) + 
                ggplot2::geom_point(size = 1.5) +
                ggplot2::scale_colour_manual(values = c(keep = "black", drop = "grey60"), limits = c("keep", "drop"), drop = FALSE, name = "Status") +
                ggplot2::theme_classic(base_size = 12) +
                ggplot2::theme(legend.position = "top") +
                ggplot2::labs(title = title, x = x_label, y = variable_name)

        # rtime use fixed y-axis and hide y-axis information.
        if (variable_name == "rtime") {
            plot = plot + 
                ggplot2::scale_y_continuous(breaks = NULL, limits = c(0.5, 1.5)) +
                ggplot2::labs(y = NULL) + 
                ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
        }

        return(plot)
    }

    plots = list()    
    ms_level_names = unique(names(record_df))

    for (ms_level in ms_level_names) {

        plots[[ms_level]] = list()

        ms_variable_names = names(record_df[[ms_level]])

        for (ms_variable in ms_variable_names) {
            data_df = record_df[[ms_level]][[ms_variable]]
            
            if (!identical(unique(as.character(data_df$variable_name)), ms_variable)) {
                stop("Unexpected variable name for ", ms_level, ".", ms_variable)
            }
            
            plot_title = paste(ms_level, ms_variable, sep = " / ")

            if (ms_variable == "rtime") {
                rt_data = data.frame(
                    spectrum_index = data_df$values,
                    values = rep(1, nrow(data_df)),
                    status = data_df$status
                )
                plots[[ms_level]][[ms_variable]] = do_plot(rt_data, variable_name = ms_variable, title = plot_title, x_label = "Retention time")
            } else {
                plot_original = do_plot(data_df, variable_name = ms_variable, title = paste(plot_title, "— spectrum order"), x_label = "Spectrum index")

                sorted_data = dplyr::arrange(data_df, values)
                sorted_data$spectrum_index = seq_len(nrow(sorted_data))

                plot_sorted = do_plot(sorted_data, variable_name = ms_variable, title = paste(plot_title, "— ascending value order"), x_label = "Sorted index")

                plots[[ms_level]][[ms_variable]] =
                    patchwork::wrap_plots(
                        plot_original,
                        plot_sorted,
                        ncol = 1,
                        guides = "collect"
                    ) & ggplot2::theme(legend.position = "top")
            }
        }
    }
    return(plots)
}
