# make_ms_qc_config <- function(
#     ms1_qc_method = 'mean_sd', 
#     ms2_qc_method = , 
#     ms1_qc_threshold,
#     ms2_qc_threshold, 
#     plot_option) {
# }



make_glycan_library_config <- function(
    glycan_type,

    mode, 

    min_charge_state,
    max_charge_state,

    monosaccharides_additional_monoisotopic_weight_vector,
    monosaccharides_additional_min_num_vector,
    monosaccharides_additional_max_num_vector,

    adduct_type,

    adduct_additional_monoisotopic_weight_vector,
    adduct_additional_charge_state_vector, 

    modification_monoisotopic_weight_vector, 
    modification_min_num_vector,
    modification_max_num_vector,

    derivatization_type,
    derivatization_additional_monoisotopic_weight,

    min_Hex_num, max_Hex_num,
    min_HexNAc_num, max_HexNAc_num,
    min_dHex_num, max_dHex_num,
    min_Neu5Ac_num, max_Neu5Ac_num,
    min_HexA_num, max_HexA_num,
    min_Neu5Gc_num, max_Neu5Gc_num,
    min_Pentose_num, max_Pentose_num,
    min_KDN_num, max_KDN_num, 
    min_EtNP_num, max_EtNP_num,
    min_AHM_num, max_AHM_num,

    min_total_monosaccharides_num,
    max_total_monosaccharides_num
) {
    structure(
        list(
            glycan_type = glycan_type,

            mode = mode, 

            min_charge_state = min_charge_state,
            max_charge_state = max_charge_state,

            custom_monosaccharide_mass = monosaccharides_additional_monoisotopic_weight_vector,
            custom_monosaccharide_min = monosaccharides_additional_min_num_vector,
            custom_monosaccharide_max = monosaccharides_additional_max_num_vector,

            adduct_type = adduct_type,

            custom_adduct_mass = adduct_additional_monoisotopic_weight_vector,
            custom_adduct_charge = adduct_additional_charge_state_vector,

            custom_modification_mass = modification_monoisotopic_weight_vector, 
            custom_modification_min = modification_min_num_vector,
            custom_modification_max = modification_max_num_vector,

            derivatization_type = derivatization_type,
            custom_derivatization_mass = derivatization_additional_monoisotopic_weight, 

            monosaccharide_range = list(
                Hex = list(min_Hex = min_Hex_num, max_Hex = max_Hex_num),
                HexNAc = list(min_HexNAc = min_HexNAc_num, max_HexNAc = max_HexNAc_num),
                dHex = list(min_dHex = min_dHex_num, max_dHex = max_dHex_num),
                Neu5Ac = list(min_Neu5Ac = min_Neu5Ac_num, max_Neu5Ac = max_Neu5Ac_num),
                HexA = list(min_HexA = min_HexA_num, max_HexA = max_HexA_num),
                Neu5Gc = list(min_Neu5Gc = min_Neu5Gc_num, max_Neu5Gc = max_Neu5Gc_num),
                Pentose = list(min_Pentose = min_Pentose_num, max_Pentose = max_Pentose_num),
                KDN = list(min_KDN = min_KDN_num, max_KDN = max_KDN_num), 
                EtNP = list(min_EtNP = min_EtNP_num, max_EtNP = max_EtNP_num),
                AHM = list(min_AHM = min_AHM_num, max_AHM = max_AHM_num)
                ),
                
                total_monosaccharide_range = list(min_total_monosaccharides = min_total_monosaccharides_num, max_total_monosaccharides = max_total_monosaccharides_num)
            ), 
          class = "glycan_library_config"
    )
}



#' Configure Adaptive MS2 Spectrum Denoising
#'
#' @description
#' Create a configuration that selects smoothing-spline or segmented
#' regression denoising independently for each MS2 spectrum, based on
#' the m/z spacing of its most intense peaks.
#' This is the default configuration used by `MS2SpectraDenoising()`.
#'
#' @inheritParams DenoisingConfigSpline
#'
#' @param top_peaks_for_pattern A single finite, positive integer-valued
#'   numeric value specifying how many of the most intense peaks are
#'   used to assess m/z spacing. Defaults to `4`.
#' @param top_peaks_evenly_spaced_tol A single finite, non-negative numeric
#'   relative tolerance for deviations of adjacent m/z gaps from their
#'   median. Defaults to `0.05`, corresponding to a relative deviation
#'   of at most 5 percent. This is not an absolute m/z or ppm tolerance.
#' @param segmented_nonlinear_fit_function A function applied to rank
#'   indices above the estimated breakpoint for an auxiliary regression
#'   fit. Defaults to `function(z) z^2 + z`. It should accept a numeric
#'   vector and return a finite numeric vector of the same length.
#'   This function does not determine the breakpoint itself.
#'   See [DenoisingConfigSegmented()].
#'
#' @details
#' For each spectrum, the most intense `top_peaks_for_pattern` peaks
#' are selected using their original intensities. Their m/z values
#' are sorted, and adjacent gaps are calculated.
#'
#' Spline denoising is selected when every gap is positive and every
#' relative deviation from the median gap is at most
#' `top_peaks_evenly_spaced_tol`. Otherwise, segmented regression is
#' selected. This tests regular spacing without requiring a particular
#' absolute gap or precursor charge.
#'
#' The spline branch uses the smoothing-parameter search, scoring
#' weights, and threshold-search settings documented in
#' [DenoisingConfigSpline()]. The segmented branch follows
#' [DenoisingConfigSegmented()]. A fitting error does not automatically
#' trigger a fallback to the other method.
#'
#' Before method selection, spectra with fewer than
#' `top_peaks_for_pattern + top_peaks_to_remove + 4` peaks are assigned
#' denoising status `"fail"` and a noise threshold of zero.
#' This check uses the original peak count, before duplicate intensity
#' values are removed for fitting. Passing it does not guarantee that
#' enough distinct intensity values remain for a successful fit.
#'
#' During peak filtering, only finite original intensities strictly
#' greater than the estimated threshold are retained. A threshold of
#' zero therefore still removes zero, negative, and non-finite
#' intensities.
#'
#' This constructor only stores settings. Configuration checks and
#' threshold estimation occur when the object is supplied to
#' `MS2SpectraDenoising()`.
#'
#' @return A list of class `"denoising_config"` containing
#'   `method = "adaptive"` and all supplied arguments as named elements.
#'
#' @seealso [DenoisingConfigSpline()], [DenoisingConfigSegmented()],
#'   `MS2SpectraDenoising()`.
#' @family MS2 denoising configuration
#' @export
#'
#' @examples
#' # Default adaptive configuration
#' config <- DenoisingConfigAdaptive()
#' config$method
#'
#' # Customize the spacing test and spline search
#' config <- DenoisingConfigAdaptive(
#'   transform = "sqrt",
#'   top_peaks_for_pattern = 5,
#'   top_peaks_evenly_spaced_tol = 0.1,
#'   spar_start = -1,
#'   spar_end = 1,
#'   spar_step = 0.05,
#'   top_peaks_to_remove = 3
#' )
DenoisingConfigAdaptive = function(
    # data transform 
    transform = "log2", 

    # determine the method for denoising
    top_peaks_for_pattern = 4, 
    top_peaks_evenly_spaced_tol = 0.05,

    # segmented
    segmented_nonlinear_fit_function = function(z) z^2+z, 
    
    # spline
    spar_start = -1.5, spar_end = 1.5, spar_step = 0.02, 
    RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1, 
    spline_cv = TRUE, 
    top_peaks_to_remove = 5
    ) {
        config = structure(
            list(
                method = "adaptive", 
                transform = transform,
                top_peaks_for_pattern = top_peaks_for_pattern, 
                top_peaks_evenly_spaced_tol = top_peaks_evenly_spaced_tol, 
                segmented_nonlinear_fit_function = segmented_nonlinear_fit_function, 
                spar_start = spar_start, spar_end = spar_end, spar_step = spar_step,
                RMSE_weight = RMSE_weight, CV_weight = CV_weight, D2_weight = D2_weight, D1_weight = D1_weight,
                spline_cv = spline_cv, 
                top_peaks_to_remove = top_peaks_to_remove
            ), 
            class = "denoising_config" 
        )
    return(config)
}


#' Configure Smoothing-Spline MS2 Spectrum Denoising
#'
#' @description
#' Create a configuration for estimating an MS2 noise threshold from
#' a smoothing spline fitted to ranked, transformed intensity values.
#' Smoothing parameters are compared using a weighted combination of
#' fit error, cross-validation, and derivative-based criteria.
#'
#' @param transform Intensity transformation used for threshold
#'   estimation. Either `"log2"` for `log2(z + 1)`, `"sqrt"` for
#'   `sqrt(z)`, `"asinh"` for `asinh(z)`, `"none"` for the identity
#'   transformation, or a custom function. Defaults to `"log2"`.
#'   A custom function should accept a numeric intensity vector and
#'   return a finite numeric vector of the same length while preserving
#'   intensity order, because transformed ranks are mapped back to
#'   sorted original intensities. Retained peak intensities are not
#'   replaced by their transformed values.
#' @param spar_start,spar_end Finite numeric scalars specifying the
#'   starting and ending values of the candidate smoothing-parameter
#'   sequence. Defaults to `-1.5` and `1.5`, respectively.
#'   Supply `spar_start <= spar_end`. Candidates are generated using
#'   `seq(spar_start, spar_end, by = spar_step)`.
#' @param spar_step A single finite, positive numeric increment between
#'   candidate smoothing parameters. Defaults to `0.02`.
#' @param RMSE_weight A single finite, non-negative numeric weight for
#'   normalized root mean squared fitting error. Defaults to `0.3`.
#' @param CV_weight A single finite, non-negative numeric weight for the
#'   normalized cross-validation criterion returned by
#'   [stats::smooth.spline()]. Defaults to `0.3`. If this criterion is
#'   unavailable for a candidate, its weight is added to `RMSE_weight`
#'   for that candidate.
#' @param D2_weight A single finite, non-negative numeric weight for
#'   normalized integrated squared second derivative, approximated on
#'   the evaluation grid. Defaults to `0.3`.
#' @param D1_weight A single finite, non-negative numeric weight for
#'   normalized variation in the first derivative, calculated as the
#'   sum of absolute differences between adjacent derivative values
#'   on the evaluation grid. Defaults to `0.1`.
#' @param spline_cv A single non-missing logical value passed to the
#'   `cv` argument of [stats::smooth.spline()]. `TRUE` requests ordinary
#'   leave-one-out cross-validation; `FALSE` requests generalized
#'   cross-validation. Defaults to `TRUE`.
#' @param top_peaks_to_remove A single finite, non-negative integer-valued
#'   numeric value specifying how many of the highest distinct-intensity
#'   ranks to exclude from the derivative-based threshold search.
#'   Defaults to `5`. These observations remain in the spline fit;
#'   this argument does not directly remove the most intense peaks
#'   from the output spectrum.
#'
#' @details
#' Duplicate original intensity values are removed before transformation.
#' The transformed values are sorted and fitted against their rank
#' indices using [stats::smooth.spline()] for each candidate `spar`.
#'
#' RMSE is evaluated at the observed ranks. Derivative criteria are
#' evaluated on a grid containing twice as many points as the number
#' of distinct original intensity values. Each scoring criterion is
#' min-max normalized across finite candidate values.
#'
#' The candidate with the smallest finite weighted score is selected;
#' ties are resolved by choosing the first candidate. Weights are
#' used directly, are not automatically normalized, and need not
#' sum to one.
#'
#' The selected spline is searched for the smallest absolute first
#' derivative at eligible grid positions. Positions above
#' `n - top_peaks_to_remove` are excluded, where `n` is the number of
#' distinct original intensity values. The selected position is rounded
#' down to an integer rank and mapped to the corresponding sorted
#' original intensity to obtain the noise threshold.
#'
#' During peak filtering, finite original intensities strictly greater
#' than the threshold are retained. Peaks equal to the threshold
#' are removed.
#'
#' Spectra with fewer than `top_peaks_to_remove + 4` original peaks
#' are assigned status `"fail"` and a threshold of zero. Duplicate
#' intensities can leave too few distinct observations for fitting
#' even when this initial check passes. Fitting or scoring errors
#' can stop denoising.
#'
#' This constructor only stores settings. Configuration checks and
#' threshold estimation occur in `MS2SpectraDenoising()`.
#'
#' @return A list of class `"denoising_config"` containing
#'   `method = "spline"` and all supplied arguments as named elements.
#'
#' @seealso [DenoisingConfigAdaptive()], [DenoisingConfigSegmented()],
#'   `MS2SpectraDenoising()`.
#' @family MS2 denoising configuration
#' @export
#'
#' @examples
#' # Default spline configuration
#' config <- DenoisingConfigSpline()
#'
#' # Use generalized cross-validation and a coarser search
#' config <- DenoisingConfigSpline(
#'   transform = "log2",
#'   spar_start = -1,
#'   spar_end = 1,
#'   spar_step = 0.05,
#'   spline_cv = FALSE,
#'   top_peaks_to_remove = 3
#' )
DenoisingConfigSpline = function(
    # data transform 
    transform = "log2", 

    # spline
    spar_start = -1.5, spar_end = 1.5, spar_step = 0.02, 
    RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1, 
    spline_cv = TRUE, 
    top_peaks_to_remove = 5
    ) {
        config = structure(
            list(
                method = "spline", 
                transform = transform, 
                spar_start = spar_start, spar_end = spar_end, spar_step = spar_step,
                RMSE_weight = RMSE_weight, CV_weight = CV_weight, D2_weight = D2_weight, D1_weight = D1_weight,
                spline_cv = spline_cv, 
                top_peaks_to_remove = top_peaks_to_remove
            ), 
            class = "denoising_config" 
        )
    return(config)
}



#' Configure Segmented-Regression MS2 Spectrum Denoising
#'
#' @description
#' Create a configuration for estimating an MS2 noise threshold from
#' a breakpoint in the relationship between ranked, transformed
#' intensity values and their rank indices.
#'
#' @inheritParams DenoisingConfigSpline
#'
#' @param segmented_nonlinear_fit_function A function applied to rank
#'   indices above the estimated breakpoint for an auxiliary regression
#'   fit. Defaults to `function(z) z^2 + z`. It should accept a numeric
#'   vector and return a finite numeric vector of the same length.
#'   It transforms the predictor in the post-breakpoint fit; it does
#'   not specify the model used to estimate the breakpoint.
#'
#' @details
#' Duplicate original intensity values are removed, the selected
#' intensity transformation is applied, and the transformed values
#' are sorted.
#'
#' An initial linear model is fitted against rank index.
#' `segmented::segmented()` then estimates a breakpoint, using one
#' third of the number of distinct intensity values as the initial
#' breakpoint location. The estimated breakpoint is rounded to an
#' integer rank and mapped to the corresponding sorted original
#' intensity to obtain the noise threshold.
#'
#' After breakpoint estimation, separate auxiliary regressions are
#' fitted below and above the breakpoint. The lower segment uses
#' untransformed rank indices. The upper segment uses rank indices
#' transformed by `segmented_nonlinear_fit_function`.
#' These fits support diagnostic output and do not redefine the
#' breakpoint-derived threshold.
#'
#' During peak filtering, only finite original intensities strictly
#' greater than the threshold are retained.
#'
#' Spectra with fewer than 10 original peaks are assigned denoising
#' status `"fail"` and a threshold of zero. This check precedes
#' removal of duplicate intensities and does not guarantee that
#' enough distinct values remain for successful regression.
#' Regression errors can stop denoising.
#'
#' This constructor only stores settings. Configuration checks and
#' threshold estimation occur in `MS2SpectraDenoising()`.
#'
#' @return A list of class `"denoising_config"` containing
#'   `method = "segmented"`, `transform`, and
#'   `segmented_nonlinear_fit_function`.
#'
#' @seealso [DenoisingConfigAdaptive()], [DenoisingConfigSpline()],
#'   `MS2SpectraDenoising()`.
#' @family MS2 denoising configuration
#' @export
#'
#' @examples
#' config <- DenoisingConfigSegmented()
#'
#' # Customize the intensity transformation and auxiliary upper fit
#' config <- DenoisingConfigSegmented(
#'   transform = "asinh",
#'   segmented_nonlinear_fit_function = function(z) z^2
#' )
DenoisingConfigSegmented = function(
    # data transform 
    transform = "log2", 

    # segmented
    segmented_nonlinear_fit_function = function(z) z^2+z
    
    ) {
        config = structure(
            list(
                method = "segmented", 
                transform = transform,
                segmented_nonlinear_fit_function = segmented_nonlinear_fit_function
            ), 
            class = "denoising_config" 
        )
    return(config)
}



#' Configure Quantile-Based MS2 Spectrum Denoising
#'
#' @description
#' Create a configuration for estimating an MS2 noise threshold from
#' a quantile of transformed, distinct intensity values within each
#' spectrum.
#'
#' @inheritParams DenoisingConfigSpline
#'
#' @param quantile_prob A single finite numeric probability in `[0, 1]`
#'   specifying the quantile used to locate the noise threshold.
#'   Defaults to `0.05`.
#'
#' @details
#' Duplicate original intensity values are removed before transformation,
#' so repeated intensities do not receive additional weight during
#' quantile estimation.
#'
#' The transformed values are sorted, and their requested quantile is
#' calculated with [stats::quantile()] using its default interpolation
#' method. The largest rank whose transformed value is less than or
#' equal to that quantile is selected. The corresponding sorted
#' original intensity becomes the noise threshold.
#'
#' The resulting threshold is therefore an observed original intensity,
#' rather than the interpolated quantile itself. With an order-preserving
#' transformation, `quantile_prob = 0` selects the minimum intensity
#' and `quantile_prob = 1` selects the maximum intensity.
#'
#' During peak filtering, only finite original intensities strictly
#' greater than the threshold are retained. Consequently, a probability
#' of zero still removes peaks at the minimum intensity, while a
#' probability of one removes all peaks. Because estimation uses
#' distinct intensity values and removes ties at the cutoff, the
#' proportion of peaks removed need not equal `quantile_prob`.
#'
#' This constructor only stores settings. The current configuration
#' validator checks that `quantile_prob` is finite and non-negative;
#' probabilities above one are rejected later by [stats::quantile()]
#' during threshold estimation.
#'
#' @return A list of class `"denoising_config"` containing
#'   `method = "quantile"`, `transform`, and `quantile_prob`.
#'
#' @seealso [DenoisingConfigFixed()], [DenoisingConfigSpline()],
#'   `MS2SpectraDenoising()`.
#' @family MS2 denoising configuration
#' @export
#'
#' @examples
#' # Default quantile configuration
#' config <- DenoisingConfigQuantile()
#'
#' # Estimate the cutoff without transforming intensities
#' config <- DenoisingConfigQuantile(
#'   transform = "none",
#'   quantile_prob = 0.1
#' )
DenoisingConfigQuantile = function(
    # data transform 
    transform = "log2", 
    
    # quantile
    quantile_prob = 0.05
    ) {
        config = structure(
            list(
                method = "quantile", 
                transform = transform, 
                quantile_prob = quantile_prob
            ), 
            class = "denoising_config" 
        )
    return(config)
}



#' Configure Fixed-Threshold MS2 Spectrum Denoising
#'
#' @description
#' Create a configuration that applies the same absolute intensity
#' threshold to every MS2 spectrum.
#'
#' @param fixed_value A single finite, non-negative numeric intensity
#'   threshold, expressed in the same units as the original peak
#'   intensities. Defaults to `30`.
#'
#' @details
#' The supplied value is used directly as the noise threshold for
#' every MS2 spectrum. No intensity transformation, normalization,
#' or spectrum-specific threshold estimation is performed.
#'
#' During peak filtering, only finite original intensities strictly
#' greater than `fixed_value` are retained. Peaks equal to the
#' threshold are removed.
#'
#' A threshold of zero retains positive finite intensities and removes
#' zero, negative, and non-finite intensities. The threshold is an
#' absolute intensity value, not a percentage of the base-peak
#' intensity or total ion current.
#'
#' This constructor only stores settings. Validation occurs when the
#' configuration is supplied to `MS2SpectraDenoising()`.
#'
#' @return A list of class `"denoising_config"` containing
#'   `method = "fixed"` and `fixed_value`.
#'
#' @seealso [DenoisingConfigQuantile()], `MS2SpectraDenoising()`.
#' @family MS2 denoising configuration
#' @export
#'
#' @examples
#' # Retain finite MS2 peak intensities strictly greater than 30
#' config <- DenoisingConfigFixed()
#'
#' # Specify a different absolute intensity threshold
#' config <- DenoisingConfigFixed(fixed_value = 100)
#'
#' # Retain only positive finite intensities
#' positive_only_config <- DenoisingConfigFixed(fixed_value = 0)
DenoisingConfigFixed = function(
    # fixed value
    fixed_value = 30
    ) {
        config = structure(
            list(
                method = "fixed", 
                fixed_value = fixed_value
            ), 
            class = "denoising_config" 
        )
    return(config)
}




#' Get Default Diagnostic Fragment m/z Values
#'
#' @description
#' Return the built-in named vector of diagnostic fragment m/z values
#' used as the default `diagnostic_frags_list` in
#' `FindMS2SpectraByDiagnosticFragments()`.
#'
#' @details
#' Vector names identify diagnostic fragments and are used as variable
#' names in the `diagnostic_frags_exp` argument of
#' `FindMS2SpectraByDiagnosticFragments()`. Names are case-sensitive.
#'
#' During spectrum screening, each fragment referenced in the expression
#' is assigned a logical presence value according to whether a denoised
#' MS2 peak matches its target m/z within the supplied mass tolerance.
#' The expression combines these presence values to determine whether
#' a spectrum is selected.
#'
#' For example, `"HexNAc & HexNAc_ProA"` requires both fragments,
#' whereas `"HexNAc | Hex"` requires at least one of the two.
#' Returning the full default vector does not require every listed
#' fragment to be present; the supplied expression determines the
#' selection conditions.
#'
#' The values are fixed target m/z values and are matched directly
#' against observed peak m/z values. They are not recalculated from
#' monosaccharide masses, glycan-library settings, precursor charges,
#' or a selected derivatization.
#'
#' This function takes no arguments and performs no spectrum screening.
#' Fragment definitions, the logical expression, and mass-error
#' tolerance are validated when passed to
#' `FindMS2SpectraByDiagnosticFragments()`.
#'
#' The returned vector can be subsetted, edited, or extended before use.
#' A custom vector must contain finite positive numeric values with
#' unique, non-empty names. Every fragment name referenced in the
#' screening expression must exist in that vector.
#'
#' @return A named numeric vector of length 12 containing the following
#'   fragment names and target m/z values, in this order:
#'
#' \tabular{lr}{
#'   Fragment name \tab Target m/z \cr
#'   `HexNAc` \tab 204.08667 \cr
#'   `HexNAc_ProA` \tab 441.2708 \cr
#'   `dHex_HexNAc_ProA` \tab 587.3287 \cr
#'   `HexNAc_HexNAc_ProA` \tab 644.3502 \cr
#'   `Hex_ProA` \tab 400.2442 \cr
#'   `Hex_Hex_ProA` \tab 562.2970 \cr
#'   `dHex_Hex_ProA` \tab 546.3021 \cr
#'   `Hex` \tab 163.06007 \cr
#'   `Bi_HexNAc` \tab 407.16607 \cr
#'   `Bisecting` \tab 1009.4824 \cr
#'   `Bisecting_dHex` \tab 1155.5403 \cr
#'   `Hex_HexNAc_ProA` \tab 603.3236
#' }
#'
#' @seealso `FindMS2SpectraByDiagnosticFragments()`,
#'   [MS2SpectraDenoising()].
#'
#' @export
#'
#' @examples
#' # Inspect the built-in fragment definitions
#' fragments <- DiagnosticFragmentsConfigDefault()
#' fragments
#'
#' # Retrieve one target m/z value
#' fragments[["HexNAc"]]
#'
#' # Select a subset for a screening rule
#' selected_fragments <- fragments[c("HexNAc", "HexNAc_ProA", "Hex")]
#' fragment_expression <- "HexNAc & (HexNAc_ProA | Hex)"
#'
#' \dontrun{
#' # denoised_data is a result from MS2SpectraDenoising().
#' # ms_data_raw provides the corresponding original MS data.
#' matches <- FindMS2SpectraByDiagnosticFragments(
#'   denoised_data = denoised_data,
#'   ms_data_raw = ms_data_raw,
#'   diagnostic_frags_list = selected_fragments,
#'   diagnostic_frags_exp = fragment_expression,
#'   mass_error_tolerance_type = "ppm_tolerance",
#'   mass_error_tolerance_val = 10
#' )
#' }
DiagnosticFragmentsConfigDefault = function() {
  c(
    HexNAc =              204.08667, 
    HexNAc_ProA =         441.2708, 
    dHex_HexNAc_ProA =    587.3287, 
    HexNAc_HexNAc_ProA =  644.3502, 
    Hex_ProA =            400.2442, 
    Hex_Hex_ProA =        562.2970, 
    dHex_Hex_ProA =       546.3021, 
    Hex =                 163.06007, 
    Bi_HexNAc =           407.16607, 
    Bisecting =           1009.4824, 
    Bisecting_dHex =      1155.5403, 
    Hex_HexNAc_ProA =     603.3236
  )
}



# MSdataPreprocessing config sup func

#' Define Filtering Conditions for One MS Level
#'
#' @description
#' Collect variable-specific filtering specifications for use as the
#' `ms1` or `ms2` argument of [MSdataPreprocessingConfig()].
#'
#' @param peaksCount Optional filtering specification for the number of
#'   peaks per spectrum. Supply [FilterMeanSD()], [FilterQuantile()], or
#'   [FilterRange()]. Defaults to `NULL`, meaning no peak-count filter.
#' @param totIonCurrent Optional filtering specification for total ion
#'   current per spectrum. Supply [FilterMeanSD()], [FilterQuantile()], or
#'   [FilterRange()]. Defaults to `NULL`, meaning no total-ion-current
#'   filter.
#' @param rtime Optional retention-time filtering specification created
#'   by [FilterRange()]. Bounds must use the same units as the input
#'   `rtime` values. Defaults to `NULL`, meaning no retention-time filter.
#'
#' @details
#' Arguments set to `NULL` are omitted from the returned list.
#' Calling this function without arguments returns an empty list,
#' which specifies no filtering conditions for the corresponding MS level.
#'
#' During preprocessing, the configured conditions are combined using
#' logical AND. Each condition is evaluated against the original data
#' for the selected MS level, rather than against the output of another
#' filter. Non-finite values fail any condition applied to that variable.
#'
#' Retention time supports only [FilterRange()]. Mean/SD and quantile
#' filtering are supported for `peaksCount` and `totIonCurrent`.
#'
#' This function collects specifications without validating them.
#' Validation occurs when the enclosing configuration is passed to
#' `MSdataPreprocessing()`.
#'
#' @return A named list containing the non-`NULL` arguments, with names
#'   drawn from `peaksCount`, `totIonCurrent`, and `rtime`.
#'
#' @seealso [MSdataPreprocessingConfig()], [FilterMeanSD()],
#'   [FilterQuantile()], [FilterRange()].
#' @family MS data preprocessing configuration
#' @export
#'
#' @examples
#' # No filtering conditions
#' MSLevelFilter()
#'
#' # Combine peak-count and total-ion-current conditions
#' ms1_filters <- MSLevelFilter(
#'   peaksCount = FilterRange(20, 1000),
#'   totIonCurrent = FilterQuantile(0.05, 1)
#' )
#'
#' config <- MSdataPreprocessingConfig(ms1 = ms1_filters)
MSLevelFilter = function(
    peaksCount = NULL,
    totIonCurrent = NULL,
    rtime = NULL
) {
    filters = list(
        peaksCount = peaksCount,
        totIonCurrent = totIonCurrent,
        rtime = rtime
        )
        
  Filter(function(x) !is.null(x), filters)
}


#' Define a Mean/SD Jump-Point Filter
#'
#' @description
#' Create a specification for filtering spectra using a lower cutoff
#' derived from gaps between sorted values. Supported for `peaksCount`
#' and `totIonCurrent` in [MSLevelFilter()].
#'
#' @param n_sd A single finite, non-negative numeric multiplier of the
#'   standard deviation of adjacent differences. May be non-integer.
#'
#' @details
#' During preprocessing, finite values of the selected variable are
#' collected separately for each MS level and sorted in ascending order.
#' The difference vector is constructed as `c(0, diff(sorted_values))`.
#'
#' A difference qualifies as a jump when it is strictly greater than
#' `mean(differences) + n_sd * stats::sd(differences)`.
#' If several jumps qualify, the largest observed value immediately
#' following a qualifying jump is selected as the cutoff.
#'
#' Spectra with finite values greater than or equal to this cutoff are
#' retained by this condition. There is no finite upper cutoff.
#' The mean and standard deviation are calculated from the difference
#' vector, not from the original observed values.
#'
#' If finite observations are available but no qualifying jump is found,
#' preprocessing stops with an error. If no finite observations are
#' available, all spectra in that MS-level subset fail this condition.
#'
#' This method is not supported for `rtime`.
#' The constructor stores `n_sd` without validating it or estimating a
#' cutoff; these operations occur in `MSdataPreprocessing()`.
#'
#' @return A list with `method = "mean_sd"` and `threshold = n_sd`.
#'
#' @seealso [MSLevelFilter()], [FilterQuantile()], [FilterRange()].
#' @family MS data preprocessing configuration
#' @export
#'
#' @examples
#' FilterMeanSD(2)
#'
#' config <- MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(peaksCount = FilterMeanSD(2)),
#'   ms2 = MSLevelFilter(totIonCurrent = FilterMeanSD(1.5))
#' )
FilterMeanSD = function(n_sd) {
  list(method = "mean_sd", threshold = n_sd)
}

#' Define a Quantile-Based Filter
#'
#' @description
#' Create a specification for retaining spectra whose variable values
#' fall between two empirical quantiles. Supported for `peaksCount`
#' and `totIonCurrent` in [MSLevelFilter()].
#'
#' @param lower A single finite numeric probability specifying the lower
#'   quantile. Must satisfy `0 <= lower < upper`.
#' @param upper A single finite numeric probability specifying the upper
#'   quantile. Must satisfy `lower < upper <= 1`.
#'
#' @details
#' During preprocessing, quantiles are calculated separately for each
#' MS level from the finite values of the selected variable in the
#' original input subset. Calculation uses [stats::quantile()] with
#' its default interpolation method (`type = 7`).
#'
#' Both resulting bounds are inclusive. Tied values and interpolated
#' quantile boundaries mean that the retained proportion need not equal
#' `upper - lower`.
#'
#' Using `lower = 0` and `upper = 1` retains all finite values for this
#' condition. Non-finite values are removed even with these probabilities.
#' If no finite observations are available, all spectra in that MS-level
#' subset fail this condition.
#'
#' This method is not supported for `rtime`.
#' The constructor stores the probabilities without validating them or
#' calculating quantiles; these operations occur in
#' `MSdataPreprocessing()`.
#'
#' @return A list with `method = "quantile"` and
#'   `threshold = c(lower, upper)`.
#'
#' @seealso [MSLevelFilter()], [FilterMeanSD()], [FilterRange()],
#'   [stats::quantile()].
#' @family MS data preprocessing configuration
#' @export
#'
#' @examples
#' # Retain values at or above the fifth percentile
#' FilterQuantile(0.05, 1)
#'
#' # Retain values between the fifth and ninety-fifth percentiles
#' config <- MSdataPreprocessingConfig(
#'   ms2 = MSLevelFilter(
#'     totIonCurrent = FilterQuantile(0.05, 0.95)
#'   )
#' )
FilterQuantile = function(lower, upper) {
  list(method = "quantile", threshold = c(lower, upper))
}

#' Define a Fixed-Range Filter
#'
#' @description
#' Create a specification for retaining spectra whose variable values
#' fall within a fixed inclusive interval. Supported for `peaksCount`,
#' `totIonCurrent`, and `rtime` in [MSLevelFilter()].
#'
#' @param start A single finite, non-negative numeric lower bound.
#'   Must be strictly less than `end`.
#' @param end A single finite, non-negative numeric upper bound.
#'   Must be strictly greater than `start`.
#'
#' @details
#' During preprocessing, a spectrum passes this condition when the
#' selected variable is finite and lies within `[start, end]`, including
#' both endpoints. Bounds are used directly and are not estimated from
#' the data. Infinite bounds are not supported.
#'
#' Bounds must use the same units as the selected variable. No unit
#' conversion is performed. For example, if `rtime` is stored in seconds,
#' express a five-minute lower bound as `5 * 60`.
#'
#' For retention-time filtering, validation additionally requires
#' `start` to be at least the minimum retention time and `end` to be at
#' most the maximum retention time in the complete input MS data.
#' These checks use the full input object, not only the MS level being
#' filtered. This observed-range restriction does not apply to
#' `peaksCount` or `totIonCurrent`.
#'
#' This constructor stores the bounds without validating them or
#' filtering spectra. Validation occurs in `MSdataPreprocessing()`.
#'
#' @return A list with `method = "start_end"` and
#'   `threshold = c(start, end)`.
#'
#' @seealso [MSLevelFilter()], [FilterMeanSD()], [FilterQuantile()].
#' @family MS data preprocessing configuration
#' @export
#'
#' @examples
#' # Specify an inclusive peak-count interval
#' FilterRange(20, 1000)
#'
#' # Specify 5 to 50 minutes for data with rtime stored in seconds.
#' # The input retention-time range must contain both endpoints.
#' config <- MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(rtime = FilterRange(5 * 60, 50 * 60)),
#'   ms2 = MSLevelFilter(rtime = FilterRange(5 * 60, 50 * 60))
#' )
FilterRange = function(start, end) {
  list(method = "start_end", threshold = c(start, end))
}


#' Create an MS Data Preprocessing Configuration
#'
#' @description
#' Combine separate MS1 and MS2 filtering specifications into a
#' configuration object for `MSdataPreprocessing()`.
#'
#' @param ms1 A list of MS1 filtering specifications created by
#'   [MSLevelFilter()]. Defaults to `MSLevelFilter()`, which specifies
#'   no filtering conditions for MS1 spectra.
#' @param ms2 A list of MS2 filtering specifications created by
#'   [MSLevelFilter()]. Defaults to `MSLevelFilter()`, which specifies
#'   no filtering conditions for MS2 spectra.
#'
#' @details
#' Filtering can be configured for peak count (`peaksCount`), total ion
#' current (`totIonCurrent`), and retention time (`rtime`).
#'
#' During preprocessing, thresholds are determined independently for MS1
#' and MS2 spectra. Within each MS level, each variable's thresholds are
#' calculated from its finite values in the original input subset,
#' before applying any other configured filters. A spectrum is retained
#' only if it passes all configured conditions for its MS level.
#'
#' Filter bounds are inclusive. Spectra with non-finite values in a
#' filtered variable are removed. Variables without a filter do not
#' affect retention. Other MS levels are preserved, and the original
#' order of retained spectra is maintained.
#'
#' This constructor stores the supplied specifications without validating
#' them or processing spectra. Configuration validation occurs when the
#' object is passed to `MSdataPreprocessing()`.
#'
#' The default configuration applies no filtering conditions, but
#' `MSdataPreprocessing()` still performs its input-data checks.
#'
#' @return A list of class `"ms_data_preprocessing_config"` with two
#'   elements, `ms1` and `ms2`, containing the supplied filtering lists.
#'
#' @seealso [MSLevelFilter()], `MSdataPreprocessing()`.
#' @family MS data preprocessing configuration
#' @export
#'
#' @examples
#' # Configure no filtering conditions
#' config <- MSdataPreprocessingConfig()
#'
#' # Configure MS1 and MS2 independently
#' config <- MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(
#'     peaksCount = FilterMeanSD(2),
#'     totIonCurrent = FilterQuantile(0.05, 1)
#'   ),
#'   ms2 = MSLevelFilter(
#'     peaksCount = FilterQuantile(0.1, 1),
#'     totIonCurrent = FilterQuantile(0.05, 1)
#'   )
#' )
#'
#' config$ms1
#' config$ms2
MSdataPreprocessingConfig = function(
    ms1 = MSLevelFilter(), 
    ms2 = MSLevelFilter()
    ) {
        config = structure(
            list(
                ms1 = ms1,
                ms2 = ms2
            ), 
            class = "ms_data_preprocessing_config" 
        )
    return(config)
}




