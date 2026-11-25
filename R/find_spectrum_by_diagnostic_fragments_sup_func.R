prepare_denoised_ms_data_ms2 = function(denoised_data, ms_data_raw) {
  if (!inherits(denoised_data, "denoised_ms_data_ms2")) {
    stop("Invalid denoised_data object. Must be of class 'denoised_ms_data_ms2'.")
  }

  denoising_info = denoised_data$denoising_result

  denoising_info = dplyr::filter(denoising_info, .data[["denoising_status"]] == "success")

  success_denoised_ms2_ids = denoising_info$ms2_spectrum_id

  # check whether success_denoised_ms2_ids in ms_data_raw
  raw_spectrum_ids = as.character(ms_data_raw[["spectrumId"]])
  if (!all(success_denoised_ms2_ids %in% raw_spectrum_ids)) {
    stop("Some successfully denoised MS2 spectra are not present in the raw MS data.")
  }

  # Filter the successfully denoised ms2 spectra 
  spectrum_ids = as.character(denoised_data$denoised_ms_data[["spectrumId"]])

  if (!all(success_denoised_ms2_ids %in% spectrum_ids)) {
    stop("Some successfully denoised MS2 spectra are missing from denoised_ms_data.")
  }

  denoised_ms_data_ms2_success = denoised_data$denoised_ms_data[spectrum_ids %in% success_denoised_ms2_ids, ]

  return(denoised_ms_data_ms2_success)
}



find_ms2_spectra_by_diagnostic_fragments = function(
      denoised_ms_data_ms2, 
      ms_data_raw,
      diagnostic_frags_list,
      diagnostic_frags_exp,
      mass_error_tolerance_type,
      mass_error_tolerance_val
    ) {
      # Extract the target diagnostic fragments based on the expression provided
      diagnostic_frags_target = diagnostic_frags_list[all.vars(parse(text = diagnostic_frags_exp))]

      # get the ms1 data
      ms1_data_raw = Spectra::filterMsLevel(ms_data_raw, 1)

      # progress bar for processing MS2 spectra      
      pb = utils::txtProgressBar(min = 0, max = length(denoised_ms_data_ms2), style = 3)
      on.exit(close(pb))

      # initialize variables
      res_list = list()
      
      ms1_rt_all = Spectra::rtime(ms1_data_raw)
      ms2_peaks_data = Spectra::peaksData(denoised_ms_data_ms2)

      get_meta <- function(name, spectra_data, missing_value) {
        n <- length(spectra_data)
        if (!name %in% Spectra::spectraVariables(spectra_data)) {
            return(rep(missing_value, n))
        }
        value <- spectra_data[[name]]
        if (length(value) != n) stop("Metadata length mismatch: ", name)
        
        return(value)
      }

      for (i in seq_along(denoised_ms_data_ms2)) {
        utils::setTxtProgressBar(pb, i)
        
        peaks_matrix = ms2_peaks_data[[i]]
        mz_values = peaks_matrix[, "mz"]
        ms2_id = denoised_ms_data_ms2$spectrumId[i]
        
        if (mass_error_tolerance_type == "ppm_tolerance") {          
          presence_map <- sapply(
            names(diagnostic_frags_target), function(name) {
              target_mz <- diagnostic_frags_target[[name]]
              any(abs((mz_values - target_mz) / target_mz) * 10^6 <= mass_error_tolerance_val)
            }
          )
        } else if (mass_error_tolerance_type == "absolute_tolerance") {
          presence_map <- sapply(
            names(diagnostic_frags_target), function(name) {
              target_mz <- diagnostic_frags_target[[name]]
              any(abs(mz_values - target_mz) <= mass_error_tolerance_val)
            }
          )
        } else {
          stop("Invalid mass_error_tolerance_type. Use 'ppm_tolerance' or 'absolute_tolerance'.")
        }

        # evaluate the diagnostic fragments expression to determine if the current spectrum matches the criteria
        is_match = base::eval(parse(text = diagnostic_frags_exp), envir = as.list(presence_map))

        # presence_list_flag <- as.list(presence_map)
        # names(presence_list_flag) <- paste0(names(presence_list_flag), "_flag")

        if (is_match) {
          current_ms2_data = denoised_ms_data_ms2[i]
          current_ms2_data_rt = current_ms2_data[['rtime']]
          candidate_index_ms1 = which(ms1_rt_all < current_ms2_data_rt)
          
          if (length(candidate_index_ms1) == 0) {
            stop("No MS1 scan with retention time earlier than MS2.", call. = FALSE)
            break
          }

          closest_index_ms1 = candidate_index_ms1[which.min(current_ms2_data_rt - ms1_rt_all[candidate_index_ms1])]
          spectra_data_target_ms1 = ms1_data_raw[closest_index_ms1]

          # Extract metadata for the closest MS1 spectrum and the current MS2 spectrum          
          new_row_spectrum_info = data.frame(
            ms2_spectrum_id = get_meta(name = 'spectrumId', spectra_data = current_ms2_data, missing_value = NA),             
            ms2_precursor_charge = get_meta(name = 'precursorCharge', spectra_data = current_ms2_data, missing_value = NA), 
            ms2_precursor_mz = get_meta(name = 'precursorMz', spectra_data = current_ms2_data, missing_value = NA),
            ms2_retention_time = get_meta(name = 'rtime', spectra_data = current_ms2_data, missing_value = NA), 
            ms1_spectrum_id = get_meta(name = 'spectrumId', spectra_data = spectra_data_target_ms1, missing_value = NA),
            ms1_retention_time = get_meta(name = 'rtime', spectra_data = spectra_data_target_ms1, missing_value = NA), 

            ms2_isolation_window_target_mz = get_meta(name = 'isolationWindowTargetMz', spectra_data = current_ms2_data, missing_value = NA),
            ms2_peaks_count = get_meta(name = 'peaksCount', spectra_data = current_ms2_data, missing_value = NA),
            ms2_total_ion_current = get_meta(name = 'totIonCurrent', spectra_data = current_ms2_data, missing_value = NA),
            ms2_base_peak_mz = get_meta(name = 'basePeakMZ', spectra_data = current_ms2_data, missing_value = NA),
            ms2_base_peak_intensity = get_meta(name = 'basePeakIntensity', spectra_data = current_ms2_data, missing_value = NA),
            ms2_collision_energy = get_meta(name = 'collisionEnergy', spectra_data = current_ms2_data, missing_value = NA),
            ms2_ionisation_energy = get_meta(name = 'ionisationEnergy', spectra_data = current_ms2_data, missing_value = NA), 
            ms2_precursor_intensity = get_meta(name = 'precursorIntensity', spectra_data = current_ms2_data, missing_value = NA)
          )
          res_list[[length(res_list) + 1]] = new_row_spectrum_info
        }
      }
      
      # check res_list
      if (length(res_list) == 0) {
        message("No MS2 spectrum was identified")
        return(NULL)
      }

      extracted_spectrum_info = dplyr::bind_rows(res_list)

      return(extracted_spectrum_info)
}


























plot_diagnostic_fragments_signal_info <- function(
    denoised_ms_data_ms2,
    ms_data_raw,
    diagnostic_frags_list,
    diagnostic_frags_exp,
    mass_error_tolerance_type,
    mass_error_tolerance_val
) {
  expr <- parse(text = diagnostic_frags_exp)
  
  # get the positive-only logical expression 
  keep_positive_conditions <- function(node, negative = FALSE) {
    if (is.symbol(node)) {
      if (negative) return(NULL)
      return(node)
    }
    
    if (!is.call(node) || !is.symbol(node[[1L]])) {
      stop(
        "Use only fragment names, &, |, ! and parentheses.",
        call. = FALSE
      )
    }
    
    op <- as.character(node[[1L]])
    
    if (op == "(") {
      return(keep_positive_conditions(node[[2L]], negative))
    }
    
    if (op == "!") {
      return(keep_positive_conditions(node[[2L]], !negative))
    }
    
    if (!op %in% c("&", "|")) {
      stop(
        "Use only fragment names, &, |, ! and parentheses.",
        call. = FALSE
      )
    }
    
    left <- keep_positive_conditions(node[[2L]], negative)
    right <- keep_positive_conditions(node[[3L]], negative)

    if (is.null(left)) return(right)
    if (is.null(right)) return(left)
    
    if (negative) {
      op <- if (op == "&") "|" else "&"
    }
    call(op, left, right)
  }
  
  positive_expr <- keep_positive_conditions(expr[[1L]])
  
  if (is.null(positive_expr)) {
    positive_expr <- quote(TRUE)
  }
  
  diagnostic_frags_exp_positive <- paste(
    deparse(positive_expr, width.cutoff = 500L),
    collapse = " "
  )
  
  # get the target spectra for positive logical expression
  positive_target_spectra <- find_ms2_spectra_by_diagnostic_fragments(
    denoised_ms_data_ms2 = denoised_ms_data_ms2,
    ms_data_raw = ms_data_raw,
    diagnostic_frags_list = diagnostic_frags_list,
    diagnostic_frags_exp = diagnostic_frags_exp_positive,
    mass_error_tolerance_type = mass_error_tolerance_type,
    mass_error_tolerance_val = mass_error_tolerance_val
  )

  # extract the target MS2 spectra for positive diagnostic conditions
  if (is.null(positive_target_spectra) || nrow(positive_target_spectra) == 0L) {
    stop("No MS2 spectra satisfy the positive diagnostic conditions.")
  }
  
  polarity <- get_logical_exp_polarity(expr[[1L]])  
  polarity <- polarity[!duplicated(names(polarity))]
  
  targets <- diagnostic_frags_list[names(polarity)]
  
  fragments <- data.frame(
    fragment = names(polarity),
    target_mz = as.numeric(unlist(targets, use.names = FALSE)),
    condition = unname(polarity),
    stringsAsFactors = FALSE
  )
  
  fragments <- fragments[order(fragments$target_mz), ]

  # extract the spectrum IDs and ensure they match the positive target spectra
  ids <- as.character(denoised_ms_data_ms2[["spectrumId"]])
  target_ids <- as.character(positive_target_spectra[["ms2_spectrum_id"]])
  all_ms2_peaks_data <- Spectra::peaksData(denoised_ms_data_ms2)
  
  if (anyNA(target_ids) || !all(target_ids %in% ids)) {
    stop(
      "Target spectrum IDs are missing or do not match the supplied spectra."
    )
  }

  plot_data_list <- vector("list", length(target_ids))

  # extract the MS2 peaks data for each target ID
  for (i in seq_along(target_ids)) {
    id = target_ids[i]
    selected <- which(ids == id)

    # Ensure exactly one MS2 spectrum is found for the current target ID
    if (length(selected) != 1L) {
      stop(
        "Expected exactly one MS2 spectrum for ID: ", id, "; found: ", length(selected)
        )
    }

    ms2_peaks_data <- all_ms2_peaks_data[[selected]]
    mz <- ms2_peaks_data[, "mz"]
    intensity <- ms2_peaks_data[, "intensity"]
    
    valid <- is.finite(mz) & is.finite(intensity)

    out = fragments
    out = dplyr::mutate(out,
                        ms2_index = selected,
                        ms2_spectrum_id = id, 
                        matched_mz = NA_real_,
                        intensity = 0,
                        mass_error = NA_real_
                        )

    for (j in seq_len(nrow(out))) {
      error <- mz - out$target_mz[j]
      
      if (mass_error_tolerance_type == "ppm_tolerance") {
        error <- error / out$target_mz[j] * 1e6
      }

      candidates <- which(
        valid & abs(error) <= mass_error_tolerance_val
      )
      
      if (length(candidates) > 0L) {
        # use the candidate with the highest intensity in the mass error tolerance window
        best <- candidates[which.max(intensity[candidates])]

        out$matched_mz[j] <- mz[best]
        out$intensity[j] <- intensity[best]
        out$mass_error[j] <- error[best]
      }
    }
    plot_data_list[[i]] <- out
  }
  
  plot_data <- dplyr::bind_rows(plot_data_list)  
  # rownames(plot_data) <- NULL

  intensity_plot_data = plot_data |>
    dplyr::filter(!is.na(.data[["matched_mz"]])) |>
    dplyr::mutate(
      intensity_minmax_norm =
        (.data[["intensity"]] - min(.data[["intensity"]])) /
        (max(.data[["intensity"]]) - min(.data[["intensity"]]))
    )

  error_plot_data = dplyr::filter(plot_data, !is.na(.data[["mass_error"]]))

  mz_breaks <- sort(unique(fragments$target_mz))
  jitter_width <- if (length(mz_breaks) > 1L) {
    min(diff(mz_breaks)) * 0.15
  } else {
    max(abs(mz_breaks) * 0.01, 0.01)
  }
  x_limits <- range(mz_breaks) + c(-jitter_width, jitter_width)

  x_scale <- function() {
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = mz_breaks,
      labels = function(x) format(x, digits = 8, trim = TRUE),
      expand = ggplot2::expansion(mult = 0.1)
    )
  }
  colour_scale <- function() {
    ggplot2::scale_colour_manual(
      values = c(positive = "black", negative = "grey50"),
      limits = c("positive", "negative"),
      labels = c("Positive_condition", "Negative_condition"),
      name = "Condition",
      drop = FALSE
    )
  }
  
  p_intensity <- ggplot2::ggplot(
    intensity_plot_data,
    ggplot2::aes(
      x = .data[["target_mz"]],
      y = .data[["intensity_minmax_norm"]]
    )
  ) +
    ggplot2::geom_line(
      ggplot2::aes(group = .data[["ms2_index"]]),
      colour = "black", linewidth = 0.4, alpha = 0.5, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data[["condition"]]),
      size = 1.5
    ) +
    x_scale() + colour_scale() +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      x = NULL,
      y = "Intensity (min-max normalized)"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  p_error <- ggplot2::ggplot(
    error_plot_data,
    ggplot2::aes(
      x = .data[["target_mz"]],
      y = .data[["mass_error"]]
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey70",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data[["condition"]]),
      size = 1.5,
      position = ggplot2::position_jitter(
        width = jitter_width,
        height = 0,
        seed = 123
      )
    ) +
    x_scale() + colour_scale() +
    ggplot2::scale_y_continuous(
      limits = c(
        -mass_error_tolerance_val,
        mass_error_tolerance_val
      ),
      breaks = seq(
        -mass_error_tolerance_val,
        mass_error_tolerance_val,
        length.out = 5L
      ),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::labs(
      x = "Diagnostic fragment m/z",
      y = mass_error_tolerance_type
    ) +
    ggplot2::theme_classic(base_size = 12)

  combined <- patchwork::wrap_plots(
    p_intensity,
    p_error,
    ncol = 1,
    heights = c(2, 1),
    guides = "collect"
  ) & ggplot2::theme(legend.position = "top")
  
  return(
    list(
      diagnostic_fragments_signal_info_plot = combined,
      diagnostic_fragments_signal_info_data = plot_data
    )
  )
}