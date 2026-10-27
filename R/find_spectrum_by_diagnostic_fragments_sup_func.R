




prepare_denoised_ms_data_ms2 = function(denoised_data, ms_data_raw) {
  if (!inherits(denoised_data, "denoised_ms_data_ms2")) {
    stop("Invalid denoised_data object. Must be of class 'denoised_ms_data_ms2'.")
  }

  denoising_info = denoised_data$denoising_result

  denoising_info = dplyr::filter(denoising_info, denoising_status == "success")

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
            ms2_tot_ion_current = get_meta(name = 'totIonCurrent', spectra_data = current_ms2_data, missing_value = NA),
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




