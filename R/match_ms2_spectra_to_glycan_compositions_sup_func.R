get_required_diagnostic_ms2_output_cols <- function() {
    c(
        "ms2_spectrum_id", 
        "ms2_precursor_charge",
        "ms2_precursor_mz",
        "ms2_retention_time",
        "ms1_spectrum_id",
        "ms1_retention_time"
    )
}


get_required_glycan_library_cols <- function(glycan_library) {
    if (!inherits(glycan_library, "glycan_library") ||
        !is.list(glycan_library)) {
        stop(
            "`glycan_library` must be an object from ConstructGlycanLibrary().",
            call. = FALSE
        )
    }

    derivatization_type <- glycan_library$derivatization_type
    
    return(
        c(
            "glycan_composition",
            "adduct_type",
            "neutral_monoisotopic_mass",
            "mz_value",
            "charge_state",
            derivatization_type
        )
    )
}



match_ms2_spectra_to_glycan_library = function(
    diagnostic_ms2_result, 
    glycan_library, 
    mass_error_tolerance_type, 
    mass_error_tolerance_val, 
    max_possible_candidates_num
  ) {
    glycan_library = glycan_library$glycan_library
    temp_top_candidates = list()

    for (i in seq_len(nrow(diagnostic_ms2_result))) {
        precursor_mz = diagnostic_ms2_result[["ms2_precursor_mz"]][i]
        precursor_charge = diagnostic_ms2_result[["ms2_precursor_charge"]][i]
        ms2_id = diagnostic_ms2_result[["ms2_spectrum_id"]][i]

        if (mass_error_tolerance_type == "ppm_tolerance") {
            all_candidates = glycan_library |> 
            dplyr::filter(
                .data[["charge_state"]] == .env$precursor_charge
            ) |> 
            dplyr::mutate(
                mass_error = 
                abs(.env$precursor_mz - .data[["mz_value"]]) / .data[["mz_value"]] * 1e6
            ) |> 
            dplyr::filter(
                .data[["mass_error"]] <= .env$mass_error_tolerance_val
            ) |> 
            dplyr::arrange(mass_error)

        } else if (mass_error_tolerance_type == "absolute_tolerance") {

            all_candidates = glycan_library |> 
            dplyr::filter(.data[["charge_state"]] == .env$precursor_charge) |> 
            dplyr::mutate(mass_error = abs(.env$precursor_mz - .data[["mz_value"]])) |> 
            dplyr::filter(.data[["mass_error"]] <= .env$mass_error_tolerance_val) |> 
            dplyr::arrange(mass_error)

        } else {
            stop("Invalid 'mass_error_tolerance_type'.")
        }

        all_candidates = dplyr::mutate(
            all_candidates, 
            ms2_spectrum_id = rep(.env$ms2_id, nrow(.env$all_candidates))
        )

        if (nrow(all_candidates) > 0) {
            top_candidates = head(all_candidates, n = max_possible_candidates_num)
            temp_top_candidates[[length(temp_top_candidates) + 1]] = top_candidates
        } else {
            message("No candidates found for MS2 spectrum ID: ", ms2_id)
        }
    }

    candidates_list = dplyr::bind_rows(temp_top_candidates)

    if (nrow(candidates_list) == 0) {
        stop("No candidates found for any MS2 spectrum.")
    }


    diagnostic_ms2_result = dplyr::select(
        diagnostic_ms2_result, 
        get_required_diagnostic_ms2_output_cols()
    )


    candidates_list = dplyr::left_join(
        candidates_list,
        diagnostic_ms2_result,
        by = "ms2_spectrum_id"
    )
    return(candidates_list)
}