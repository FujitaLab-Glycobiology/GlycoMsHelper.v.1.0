validate_glycan_library_config <- function(config) {

    # check the class
    if (!inherits(config, "glycan_library_config")) {
        stop(
            "`config` must be created by make_glycan_library_config().", 
        )
    }

    # check the list structure
    if (!is.list(config)) {
        stop("The input config must be a list.")
    }

    # # check the required fields
    required_fields = glycan_lib_fields()
    for (field in required_fields) {
        if (!exists(field, where = config)) {
            stop(paste("The required field `", field, "` is missing.", sep = ""))
        }
    }

    # check the glycan type
    supported_glycan_type = get_supported_glycan_type()
    if (!is.character(config$glycan_type)) {
        stop(paste(
            "The `glycan_type` must be a character string.\n", 
            paste0("Build in `glycan_type`: ", paste0("`", supported_glycan_type, "`", collapse = ", "))))
    }
    if (!config$glycan_type %in% supported_glycan_type) {
        stop(paste(
            "The `glycan_type` is not supported.\n", 
            paste0("Build in `glycan_type`: ", paste0("`", supported_glycan_type, "`", collapse = ", "))))
    }


    # check the mode
    if (!config$mode %in% c("positive", "negative")) {
        stop(paste(
            "The `mode` must be a character string.\n", 
            "Build in `mode`: `positive`, `negative`"))
    }

    # check the charge state
    # positive integer
    validate_single_positive_number(config$min_charge_state, "min_charge_state", is_integers = TRUE)
    validate_single_positive_number(config$max_charge_state, "max_charge_state", is_integers = TRUE)

    # max >= min
    if (config$min_charge_state > config$max_charge_state) {
        stop("The `max_charge_state` must be greater than or equal to the `min_charge_state`.")
    }



    # initialize flags
    config$flags = list(
        has_custom_monosaccharide = FALSE,
        has_custom_adduct = FALSE,
        has_custom_modification = FALSE,
        has_custom_derivatization = FALSE
    )

    # check the custom monosaccharides
    has_custom_monosaccharide = !is.null(config$custom_monosaccharide_mass) || 
        !is.null(config$custom_monosaccharide_min) ||
        !is.null(config$custom_monosaccharide_max)

    complete_custom_monosaccharide = !is.null(config$custom_monosaccharide_mass) &&
        !is.null(config$custom_monosaccharide_min) &&
        !is.null(config$custom_monosaccharide_max)

    # biult-in monosaccharides 
    built_in_monosaccharides = get_reference_monosaccharides()

    if (complete_custom_monosaccharide) {
        # not allow built-in monosaccharides in customized monosaccharides
        if(any(names(config$custom_monosaccharide_mass) %in% names(built_in_monosaccharides))) {
            stop(paste0(
                "The `monosaccharides_additional_monoisotopic_weight_vector` contains name conflicts with built-in monosaccharides.\n",
                "Please provide only customized monosaccharides.")
            )
        }

        # positive numeric
        validate_positive_number_vector(
            x = config$custom_monosaccharide_mass, 
            name = "monosaccharides_additional_monoisotopic_weight_vector", 
            is_integers = FALSE, 
            is_named_vector = TRUE)

        # full names
        if(
            !setequal(names(config$custom_monosaccharide_mass), names(config$custom_monosaccharide_min)) ||
            !setequal(names(config$custom_monosaccharide_mass), names(config$custom_monosaccharide_max)) ||
            !(length(config$custom_monosaccharide_mass) == length(config$custom_monosaccharide_min)) ||
            !(length(config$custom_monosaccharide_mass) == length(config$custom_monosaccharide_max))
        ) {
            stop(paste(
                "The\n", 
                "`monosaccharides_additional_monoisotopic_weight_vector`\n", 
                "`monosaccharides_additional_min_num_vector`\n", 
                "`monosaccharides_additional_max_num_vector`\n", 
                "must have the same names and lengths."
            ))
        } 

        # non negative integers named vector
        validate_nonnegative_number_vector(config$custom_monosaccharide_min, "monosaccharides_additional_min_num_vector", is_integers = TRUE, is_named_vector = TRUE)
        validate_nonnegative_number_vector(config$custom_monosaccharide_max, "monosaccharides_additional_max_num_vector", is_integers = TRUE, is_named_vector = TRUE)

        # max >= min
        for (cust_mono in names(config$custom_monosaccharide_min)) {
            if (config$custom_monosaccharide_min[cust_mono] > config$custom_monosaccharide_max[cust_mono]) {
                stop(paste("The `monosaccharides_additional_max_num_vector` for", cust_mono, "must be greater than or equal to the `monosaccharides_additional_min_num_vector`."))
            }
        }

        # can not define like this: monosaccharides_additional_min_num_vector = c(example = c(1, 2))
        for (cust_mono in names(config$custom_monosaccharide_min)) {
            validate_single_nonnegative_number(
                x = config$custom_monosaccharide_min[[cust_mono]],
                name = paste0("min_", cust_mono),
                is_integers = TRUE
            )
        }
        for (cust_mono in names(config$custom_monosaccharide_max)) {
            validate_single_nonnegative_number(
                x = config$custom_monosaccharide_max[[cust_mono]],
                name = paste0("max_", cust_mono),
                is_integers = TRUE
            )
        }

        # flag
        config$flags$has_custom_monosaccharide = TRUE

    } else if (has_custom_monosaccharide && !complete_custom_monosaccharide) {
            stop(paste(
                "If you want to customize the monosaccharides, you must provide all three parameters:\n", 
                "`monosaccharides_additional_monoisotopic_weight_vector`\n", 
                "`monosaccharides_additional_min_num_vector`\n", 
                "`monosaccharides_additional_max_num_vector`"
            ))
    }





    # check the adduct
    has_customized_adduct = !is.null(config$custom_adduct_mass) || !is.null(config$custom_adduct_charge)    
    complete_customized_adduct = !is.null(config$custom_adduct_mass) && !is.null(config$custom_adduct_charge)
    
    # built-in adducts
    built_in_adducts = get_reference_adducts()

    # customized adducts
    if (complete_customized_adduct) {
        # not allow built-in adducts in customized adducts
        if (any(names(config$custom_adduct_mass) %in% built_in_adducts$adduct)) {
            stop(paste0(
                "The `adduct_additional_monoisotopic_weight_vector` contains built-in adducts. \n",
                "Please provide only customized adducts.")
            )
        }

        # positive numeric
        validate_positive_number_vector(
            x = config$custom_adduct_mass, 
            name = "adduct_additional_monoisotopic_weight_vector", 
            is_integers = FALSE, 
            is_named_vector = TRUE)

        # full names
        if(
            !setequal(names(config$custom_adduct_mass), names(config$custom_adduct_charge)) ||
            !(length(config$custom_adduct_mass) == length(config$custom_adduct_charge))
        ) {
            stop(paste(
                "The `adduct_additional_monoisotopic_weight_vector` and `adduct_additional_charge_state_vector` must have the same names and lengths."
            ))
        } 

        # positive integers named vector
        validate_positive_number_vector(
            x = config$custom_adduct_charge, 
            name = "adduct_additional_charge_state_vector", 
            is_integers = TRUE, 
            is_named_vector = TRUE)

        # not allow built-in adducts in customized adducts
        if (any(names(config$custom_adduct_mass) %in% built_in_adducts$adduct)) {
            stop(paste0(
                "The `adduct_additional_monoisotopic_weight_vector` contains built-in adducts. ",
                "Please provide only customized adducts.")
            )
        }

        # duplicate adducts check
        if (any(duplicated(names(config$custom_adduct_mass))) || any(duplicated(names(config$custom_adduct_charge)))) {
            stop("The `adduct_additional_monoisotopic_weight_vector` or `adduct_additional_charge_state_vector` contains duplicate entries.")
        }

        # flag
        config$flags$has_custom_adduct = TRUE

    } else if (has_customized_adduct && !complete_customized_adduct) {
            stop(paste(
                "If you want to customize the adducts, you must provide both parameters:\n", 
                "`adduct_additional_monoisotopic_weight_vector` and `adduct_additional_charge_state_vector`."))
    }

    # build in adducts validation
    if ((is.null(config$adduct_type) || length(config$adduct_type) == 0) && !complete_customized_adduct) {
        stop("The `adduct_type` must not be empty with no customized adducts.")
    } else if (length(config$adduct_type) == 0 && !complete_customized_adduct) {
        stop("The `adduct_type` must not be empty with no customized adducts.")
    } else if (length(config$adduct_type) == 0 && complete_customized_adduct) {
        # allow empty adduct_type if customized adducts are complete

        # calculate the minimum charge state of the adducts 
        adduct_charge_state_min = min(config$custom_adduct_charge)

    } else if (length(config$adduct_type) >= 1) {

        if (!all(config$adduct_type %in% built_in_adducts$adduct)) {
            stop(paste0(
                "The `adduct_type` must be in the following:\n",
                "  - Built-in adducts: ",
                paste(built_in_adducts$adduct, collapse = ", "),
                "\n",
                "  - Customized adducts should provided in `adduct_additional_monoisotopic_weight_vector` ",
                "and `adduct_additional_charge_state_vector`."
                )
            )
        }

        # duplicate adducts check
        if (any(duplicated(config$adduct_type))) {
            stop("The `adduct_type` contains duplicate entries.")
        }

        # calculate the minimum charge state of the adducts 
        built_in_adduct_charge_state_min = min(sapply(config$adduct_type, function(x) built_in_adducts$charge[built_in_adducts$adduct == x]))

        if (complete_customized_adduct) {
            customized_adduct_charge_state_min = min(config$custom_adduct_charge)
            adduct_charge_state_min = min(built_in_adduct_charge_state_min, customized_adduct_charge_state_min)
        } else {
            adduct_charge_state_min = built_in_adduct_charge_state_min
        }
    }

    # max_charge_state >= minimum charge state of the adducts
    if (adduct_charge_state_min > config$max_charge_state) {
        stop(paste("The `max_charge_state` must be greater than or equal to the minimum charge state of the adducts."))
    }

    # customized modification  
    has_customized_modification = !is.null(config$custom_modification_mass) || 
    !is.null(config$custom_modification_min) ||
    !is.null(config$custom_modification_max)
    
    complete_customized_modification = !is.null(config$custom_modification_mass) &&
    !is.null(config$custom_modification_min) &&
    !is.null(config$custom_modification_max)

    # check
    if (complete_customized_modification) {
        # positive numeric
        validate_positive_number_vector(
            x = config$custom_modification_mass, 
            name = "modification_monoisotopic_weight_vector", 
            is_integers = FALSE, 
            is_named_vector = TRUE)

        # full names
        if(
            !setequal(names(config$custom_modification_mass), names(config$custom_modification_min)) ||
            !setequal(names(config$custom_modification_mass), names(config$custom_modification_max)) ||
            !(length(config$custom_modification_mass) == length(config$custom_modification_min)) ||
            !(length(config$custom_modification_mass) == length(config$custom_modification_max))
        ) {
            stop("The `modification_monoisotopic_weight_vector`, `modification_min_num_vector`, and `modification_max_num_vector` must have the same names and lengths.")
        } 

        # non negative integers named vector
        validate_nonnegative_number_vector(config$custom_modification_min, "modification_min_num_vector", is_integers = TRUE, is_named_vector = TRUE)
        validate_nonnegative_number_vector(config$custom_modification_max, "modification_max_num_vector", is_integers = TRUE, is_named_vector = TRUE)

        # max >= min
        for (modi in names(config$custom_modification_min)) {
            if (config$custom_modification_min[modi] > config$custom_modification_max[modi]) {
                stop(paste("The `modification_max_num_vector` for", modi, "must be greater than or equal to the `modification_min_num_vector`."))
            }
        }

        # check for duplicate modification names
        if (any(duplicated(names(config$custom_modification_mass))) ||
            any(duplicated(names(config$custom_modification_min))) ||
            any(duplicated(names(config$custom_modification_max)))) {
            stop("The `modification_monoisotopic_weight_vector`, `modification_min_num_vector`, or `modification_max_num_vector` contains duplicate entries.")
        }

        # flag
        config$flags$has_custom_modification = TRUE

    } else if (has_customized_modification && !complete_customized_modification) {
            stop(paste(
                "If you want to customize the modifications, you must provide all three parameters:\n", 
                "`modification_monoisotopic_weight_vector`\n", 
                "`modification_min_num_vector`\n", 
                "`modification_max_num_vector`"))
    }

    # check the derivatization
    built_in_derivatizations = get_reference_derivatizations()

    # Both built-in and custom derivatizations require a valid name
    if (!is.character(config$derivatization_type) || 
        length(config$derivatization_type) != 1L || 
        is.na(config$derivatization_type) || 
        !nzchar(trimws(config$derivatization_type))) {
            stop(
                paste0(
                "`derivatization_type` must be a single non-empty string.\n",
                "For a built-in derivatization, provide one of: ",
                paste(names(built_in_derivatizations), collapse = ", "),
                ".\n",
                "For a custom derivatization, provide a new name and ",
                "`derivatization_additional_monoisotopic_weight`."
                )
            )
    }







    if (config$derivatization_type %in% names(built_in_derivatizations)) {
         # Built-in derivatization: use its reference mass
         if (!is.null(config$custom_derivatization_mass)) {
            stop(
                paste0(
                    "Do not provide ",
                    "`derivatization_additional_monoisotopic_weight` ",
                    "when using a built-in derivatization.\n",
                    "To use a custom mass, provide a new ",
                    "`derivatization_type` name."
                )
            )
        }
    } else {
        # Custom derivatization: require a mass
        if (is.null(config$custom_derivatization_mass)) {
            stop(
                paste0(
                    "A custom `derivatization_type` requires ",
                    "`derivatization_additional_monoisotopic_weight`."
                )
            )
        }
        
        validate_single_nonnegative_number(
            x = config$custom_derivatization_mass, 
            name = "derivatization_additional_monoisotopic_weight", 
            is_integers = FALSE
        )

        # flag
        config$flags$has_custom_derivatization = TRUE
    }



    validate_count_range = function(count_range, min_name, max_name) {
        if (length(count_range) != 2L) {
            stop(
                paste0(
                    "`", min_name, "` and `", max_name,
                    "` must each be supplied."
                ),
                call. = FALSE
            )
        }

        if (is.list(count_range)) {
            expected_names = sub("_num$", "", c(min_name, max_name))

            if (!identical(names(count_range), expected_names)) {
                stop(
                    paste0(
                        "A count range must contain named `", expected_names[1], "` and `", expected_names[2], "` entries."
                    ),
                    call. = FALSE
                )
            }
        } else if (!is.numeric(count_range)) {
            stop(
                "A normalized count range must be numeric.",
                call. = FALSE
            )
        }

        min_value = count_range[[1]]
        max_value = count_range[[2]]

        validate_single_nonnegative_number(x = min_value, name = min_name, is_integers = TRUE)

        validate_single_nonnegative_number(x = max_value, name = max_name, is_integers = TRUE)

        if (min_value > max_value) {
            stop(
                paste0(
                    "`", max_name,
                    "` must be greater than or equal to `",
                    min_name, "`."
                ),
                call. = FALSE
            )
        }
        return(c(min_value, max_value))
    }

    # Check each built-in monosaccharide range.
    for (mono_name in names(config$monosaccharide_range)) {
        config$monosaccharide_range[[mono_name]] = validate_count_range(
            count_range = config$monosaccharide_range[[mono_name]],
            min_name = paste0("min_", mono_name),
            max_name = paste0("max_", mono_name)
        )
    }

    # Check total monosaccharide range.
    config$total_monosaccharide_range = validate_count_range(
        count_range = config$total_monosaccharide_range,
        min_name = "min_total_monosaccharides",
        max_name = "max_total_monosaccharides"
    )



    # check the total monosaccharide range against the sum of individual monosaccharide ranges
    built_in_mono_total_min_sum = sum(sapply(config$monosaccharide_range, function(x) x[1]))
    built_in_mono_total_max_sum = sum(sapply(config$monosaccharide_range, function(x) x[2]))

    # sum
    if (complete_custom_monosaccharide) {
        mono_total_min_sum = built_in_mono_total_min_sum + sum(config$custom_monosaccharide_min)
        mono_total_max_sum = built_in_mono_total_max_sum + sum(config$custom_monosaccharide_max)
    } else {
        mono_total_min_sum = built_in_mono_total_min_sum
        mono_total_max_sum = built_in_mono_total_max_sum
    }

    # max total monosaccharide range must be greater than or equal to the sum of the min values for individual monosaccharides
    if (mono_total_min_sum > config$total_monosaccharide_range[2]) {
        stop(paste("The `max_total_monosaccharides_num` for total monosaccharides must be greater than or equal to the sum of the `min` values for individual monosaccharides."))        
    }

    # min total monosaccharide range must be less than or equal to the sum of the max values for individual monosaccharides
    if (mono_total_max_sum < config$total_monosaccharide_range[1]) {
        stop(paste("The `min_total_monosaccharides_num` for total monosaccharides must be less than or equal to the sum of the `max` values for individual monosaccharides."))
    }

    validate_single_nonnegative_number(
        x = config$total_monosaccharide_range[1],
        name = "min_total_monosaccharides_num",
        is_integers = TRUE
    )
    validate_single_nonnegative_number(
        x = config$total_monosaccharide_range[2],
        name = "max_total_monosaccharides_num",
        is_integers = TRUE
    )



    # check the validity of column names for built-in and custom components, modifications, derivatization, and derived columns
    name_groups <- list(
        built_in_monosaccharides = names(get_reference_monosaccharides()),
        custom_monosaccharides = names(config$custom_monosaccharide_mass),

        built_in_adducts = get_reference_adducts()$adduct,
        custom_adducts = names(config$custom_adduct_mass),

        modifications = names(config$custom_modification_mass),
        derivatization = config$derivatization_type,

        derived_columns = c(
            "total_monosaccharides",
            "charge_state",
            "neutral_monoisotopic_mass",
            "mz_value"
        )
    )

    for (group_name in names(name_groups)) {
        group_names <- name_groups[[group_name]]

        if (anyNA(group_names) ||
            any(!nzchar(trimws(group_names)))) {
            stop(
                paste0("Invalid column names in `", group_name, "`."),
                call. = FALSE
            )
        }
    }

    all_names <- unlist(name_groups, use.names = FALSE)
    conflicts <- unique(all_names[duplicated(all_names)])

    if (length(conflicts) > 0L) {
        details <- vapply(
            conflicts,
            function(column_name) {
                owners <- names(name_groups)[
                    vapply(
                        name_groups,
                        function(x) column_name %in% x,
                        logical(1)
                    )
                ]

                paste0(
                    "  - `", column_name, "`: ",
                    paste(owners, collapse = ", ")
                )
            },
            character(1)
        )

        stop(
            paste(
                c(
                    "Column name conflicts detected:",
                    details,
                    "Please rename the customized components."
                ),
                collapse = "\n"
            ),
            call. = FALSE
        )
    }

    return(config)

}
   













































validate_denoising_config <- function(config) {

    # check the class
    if (!inherits(config, "denoising_config")) {
        stop(
            paste(
                "`config` must be created by one of the denoising config functions: \n", 
                "  - DenoisingConfigAdaptive\n",
                "  - DenoisingConfigSpline\n",
                "  - DenoisingConfigSegmented\n",
                "  - DenoisingConfigQuantile\n",
                "  - DenoisingConfigFixed"
            )
        )
    }

    # check the method for denoising
    built_in_denoising_methods = get_denoising_methods()

    if (!(config$method %in% built_in_denoising_methods)) {
        stop(paste(
            "The `method` must be a character string.\n", 
            paste(built_in_denoising_methods, collapse = ", ")
            )
        )
    }
    
    # determine the method for denoising
    built_in_transform_methods = get_transform_methods()

    # adaptive method
    if (config$method == "adaptive") {
        # check the transform
        if (is.character(config$transform)) {
            if (!(config$transform %in% built_in_transform_methods)) {
                stop(paste(
                    "The `transform` must be a character string.\n", 
                    paste(built_in_transform_methods, collapse = ", ")
                    )
                )
            }
        } else if (is.function(config$transform)) {
            # do nothing, user-defined function
        } else {
            stop(paste(
                "The `transform` must be a character string or a function.\n", 
                "Build in `transform`: ", paste(built_in_transform_methods, collapse = ", ")
                )
            )
        }

        # check the top_peaks_for_pattern
        validate_single_positive_number(config$top_peaks_for_pattern, "top_peaks_for_pattern", is_integers = TRUE)

        # check the top_peaks_evenly_spaced_tol
        validate_single_nonnegative_number(config$top_peaks_evenly_spaced_tol, "top_peaks_evenly_spaced_tol", is_integers = FALSE)
        
        # check the segmented_nonlinear_fit_function
        if (!(is.function(config$segmented_nonlinear_fit_function))) {
            stop(paste(
                "Please define a non-linear fitting function for the segmented regression"
                )
            )
        }

        # check the spline parameters
        validate_single_number(config$spar_start, "spar_start", is_integers = FALSE) 
        validate_single_number(config$spar_end, "spar_end", is_integers = FALSE)
        validate_single_positive_number(config$spar_step, "spar_step", is_integers = FALSE)
        validate_single_nonnegative_number(config$RMSE_weight, "RMSE_weight", is_integers = FALSE)
        validate_single_nonnegative_number(config$CV_weight, "CV_weight", is_integers = FALSE)
        validate_single_nonnegative_number(config$D2_weight, "D2_weight", is_integers = FALSE)
        validate_single_nonnegative_number(config$D1_weight, "D1_weight", is_integers = FALSE)

        # check the spline_cv
        validate_single_logical(config$spline_cv, "spline_cv")

        # check the top_peaks_to_remove
        validate_single_nonnegative_number(config$top_peaks_to_remove, "top_peaks_to_remove", is_integers = TRUE)

    # spline method
    } else if (config$method == "spline") {
        # check the transform
        if (is.character(config$transform)) {
            if (!(config$transform %in% built_in_transform_methods)) {
                stop(paste(
                    "The `transform` must be a character string.\n", 
                    paste(built_in_transform_methods, collapse = ", ")
                    )
                )
            }
        } else if (is.function(config$transform)) {
            # do nothing, user-defined function
        } else {
            stop(paste(
                "The `transform` must be a character string or a function.\n", 
                "Build in `transform`: ", paste(built_in_transform_methods, collapse = ", ")
                )
            )
        }

        # check the spline parameters
        validate_single_number(config$spar_start, "spar_start", is_integers = FALSE) 
        validate_single_number(config$spar_end, "spar_end", is_integers = FALSE)
        validate_single_positive_number(config$spar_step, "spar_step", is_integers = FALSE)
        validate_single_nonnegative_number(config$RMSE_weight, "RMSE_weight", is_integers = FALSE)
        validate_single_nonnegative_number(config$CV_weight, "CV_weight", is_integers = FALSE)
        validate_single_nonnegative_number(config$D2_weight, "D2_weight", is_integers = FALSE)
        validate_single_nonnegative_number(config$D1_weight, "D1_weight", is_integers = FALSE)

        # check the spline_cv
        validate_single_logical(config$spline_cv, "spline_cv")

        # check the top_peaks_to_remove
        validate_single_nonnegative_number(config$top_peaks_to_remove, "top_peaks_to_remove", is_integers = TRUE)

    # segmented method
    } else if (config$method == "segmented") {
        # check the transform
        if (is.character(config$transform)) {
            if (!(config$transform %in% built_in_transform_methods)) {
                stop(paste(
                    "The `transform` must be a character string.\n", 
                    paste(built_in_transform_methods, collapse = ", ")
                    )
                )
            }
        } else if (is.function(config$transform)) {
            # do nothing, user-defined function
        } else {
            stop(paste(
                "The `transform` must be a character string or a function.\n", 
                "Build in `transform`: ", paste(built_in_transform_methods, collapse = ", ")
                )
            )
        }

        # check the segmented_nonlinear_fit_function
        if (!(is.function(config$segmented_nonlinear_fit_function))) {
            stop(paste(
                "Please define a non-linear fitting function for the segmented regression"
                )
            )
        }

    # quantile method
    } else if (config$method == "quantile") {
        # check the quantile
        validate_single_nonnegative_number(config$quantile_prob, "quantile_prob", is_integers = FALSE)

    # fixed method
    } else if (config$method == "fixed") {
        # check the fixed value
        validate_single_nonnegative_number(config$fixed_value, "fixed_value", is_integers = FALSE)
    }

}








validate_plot_ms2_denoising_info = function(ms2_spectrum_id, ms_data, denoised_data) {
    ms2_data = Spectra::filterMsLevel(ms_data, msLevel = 2)
    ms2_spectrum_ids = as.character(ms2_data[["spectrumId"]])

    # check ms2_spectrum_id
    if (!is.character(ms2_spectrum_id) || length(ms2_spectrum_id) != 1L) {
        stop("`ms2_spectrum_id` must be a single character string (example: `function=2 process=0 scan=2`).")
    }

    if (!(ms2_spectrum_id %in% ms2_spectrum_ids)) {
        stop("The specified `ms2_spectrum_id` does not exist in the provided `ms_data`.")
    }

    # check if corresponding ms2 spectrum been denoised successfully
    denoising_status = dplyr::filter(denoised_data$denoising_result, .data$ms2_spectrum_id == .env$ms2_spectrum_id)$denoising_status
    if (length(denoising_status) == 0L) {
        stop("The specified `ms2_spectrum_id` does not have a corresponding entry in the `denoising_result`.")
    }
    if (denoising_status != "success") {
        stop("The specified `ms2_spectrum_id` has not been successfully denoised.")
    }

    # check class 
    if (!inherits(denoised_data, "denoised_ms_data_ms2")) {
        stop("`denoised_data` must be of class 'denoised_ms_data_ms2' from MS2SpectrumDenoising().")
    }

    # check the denoised_data
    if (!is.list(denoised_data) || !all(c("denoising_result", "denoised_ms_data", "config") %in% names(denoised_data))) {
        stop("`denoised_data` must be a list containing 'denoising_result', 'denoised_ms_data', and 'config' from MS2SpectrumDenoising().")
    }

    if (nrow(denoised_data$denoising_result) <= 0) {
        stop("The `denoising_result` in `denoised_data` must contain at least one row.")
    }

    # check if denoised_data$denoising_result matches the ms2 id in ms_data
    if (!identical(ms2_spectrum_ids, denoised_data$denoising_result$ms2_spectrum_id)) {
        stop("The specified `ms2_spectrum_id` does not exist in the `denoising_result` of the provided `denoised_data`.")
    }

    # check if the specified ms2_spectrum_id exists in both ms_data and denoised_data
    if (!(ms2_spectrum_id %in% ms2_spectrum_ids)) {
        stop("The specified `ms2_spectrum_id` does not exist in the provided `ms_data`.")
    }

    if (!(ms2_spectrum_id %in% denoised_data$denoising_result$ms2_spectrum_id)) {
        stop("The specified `ms2_spectrum_id` does not exist in the `denoising_result` of the provided `denoised_data`.")
    }
}


















validate_ms_data_preprocessing_config = function(ms_data, config) {
    if (!inherits(config, "ms_data_preprocessing_config")) {
        stop("`config` must be of class 'ms_data_preprocessing_config' from MSdataPreprocessingConfig().")
    }

    # validate threshold for 'mean_sd'
    validate_mean_sd_threshold = function(threshold) {
        validate_single_nonnegative_number(
            x = threshold, 
            name = "mean_sd threshold", 
            is_integers = FALSE
        )
    }
    # validate threshold for 'quantile'
    validate_quantile_threshold = function(threshold) {
        if (length(threshold) != 2) {
            stop("The 'quantile threshold' must be a vector of length 2.")
        }

        validate_nonnegative_number_vector(
            x = threshold, 
            name = "quantile threshold", 
            is_integers = FALSE, 
            is_named_vector = FALSE
        )

        if (any(threshold < 0) || any(threshold > 1)) {
            stop("All elements of 'quantile threshold' must be between 0 and 1.")
        }

        if (threshold[1] >= threshold[2]) {
            stop("The first element of 'quantile threshold' must be less than the second element.")
        }
    }
    # validate threshold for 'start_end'
    validate_start_end_threshold = function(threshold, ms_data, is_rt) {
        if (length(threshold) != 2) {
            stop("The 'start_end threshold' must be a vector of length 2.")
        }

        validate_nonnegative_number_vector(
            x = threshold, 
            name = "start_end threshold", 
            is_integers = FALSE, 
            is_named_vector = FALSE
        )

        if (threshold[1] >= threshold[2]) {
            stop("The first element of 'start_end threshold' must be less than the second element.")
        }

        if (is_rt) {
            rt_ms_data = ms_data[["rtime"]]
            if (threshold[2] > max(rt_ms_data)) {
                stop("The second element of 'start_end threshold' for retention time must not exceed the maximum retention time in the MS data.")
            }
            if (threshold[1] < min(rt_ms_data)) {
                stop("The first element of 'start_end threshold' for retention time must not be less than the minimum retention time in the MS data.")
            }
        }
    }

    # strart validating the configuration for each MS level and variable
    supported_ms_data_filtering_variables = get_supported_ms_data_filtering_variables()
    supported_ms_data_filtering_methods = get_ms_data_filtering_methods()

    # get the names of ms level in the config
    if (length(config) != 2L || !setequal(names(config), c("ms1", "ms2"))) {
        stop("config must contain exactly 'ms1' and 'ms2'.")
    }

    ms_levels = names(config)

    for (ms_level in ms_levels) {
        if (!is.list(config[[ms_level]])) {
            stop(paste0("The configuration for MS level ", ms_level, " must be a list."))
        }
        # if (length(config[[ms_level]]) == 0L || is.null(config[[ms_level]])) {
        #     stop(paste0("The configuration for MS level ", ms_level, " must not be empty."))
        # }
    }    

    for (ms_level in ms_levels) {
        ms_variables_names = names(config[[ms_level]])

        if (length(config[[ms_level]]) > 0L && 
                (is.null(ms_variables_names) || 
                anyNA(ms_variables_names) || 
                any(ms_variables_names == "") || 
                anyDuplicated(ms_variables_names) > 0L)) {
            stop("Filtering variables must have non-empty, unique names.")
        }

        if (!all(ms_variables_names %in% supported_ms_data_filtering_variables)) {
            stop(paste0("The configuration for MS level ", ms_level, " contains unsupported variables."))
        }

        for (ms_variable_name in ms_variables_names) {
            if (!is.list(config[[ms_level]][[ms_variable_name]])) {
                stop(paste0("The configuration for MS level ", ms_level, " and variable ", ms_variable_name, " must be a list."))
            }

            method_info = config[[ms_level]][[ms_variable_name]]

            if (!all(names(method_info) %in% c("method", "threshold"))) {
                stop(paste0("The configuration for MS level ", ms_level, " and variable ", ms_variable_name, " must contain only 'method' and 'threshold' fields."))
            }
            if (!("method" %in% names(method_info))) {
                stop(paste0("The configuration for MS level ", ms_level, " and variable ", ms_variable_name, " must contain a 'method' field."))
            }
            if (!("threshold" %in% names(method_info))) {
                stop(paste0("The configuration for MS level ", ms_level, " and variable ", ms_variable_name, " must contain a 'threshold' field."))
            }

            if (!is.character(method_info[["method"]]) ||
                length(method_info[["method"]]) != 1L ||
                is.na(method_info[["method"]]) ||
                !(method_info[["method"]] %in% supported_ms_data_filtering_methods)) {
                stop(paste0("The method '", method_info[["method"]], "' for MS level ", ms_level, " and variable ", ms_variable_name, " is not supported."))
            }

            if (method_info[["method"]] == "mean_sd" && ms_variable_name != "rtime") {
                validate_mean_sd_threshold(threshold = method_info[["threshold"]])
            } else if (method_info[["method"]] == "mean_sd" && ms_variable_name == "rtime") {
                stop(paste0("The 'mean_sd' method is not supported for the 'rtime' variable at MS level ", ms_level, "."))

            } else if (method_info[["method"]] == "quantile" && ms_variable_name != "rtime") {
                validate_quantile_threshold(threshold = method_info[["threshold"]])
            } else if (method_info[["method"]] == "quantile" && ms_variable_name == "rtime") {
                stop(paste0("The 'quantile' method is not supported for the 'rtime' variable at MS level ", ms_level, "."))

            } else if (method_info[["method"]] == "start_end" && ms_variable_name == "rtime") {
                validate_start_end_threshold(threshold = method_info[["threshold"]], ms_data = ms_data, is_rt = TRUE)
            } else if (method_info[["method"]] == "start_end" && ms_variable_name != "rtime") {
                validate_start_end_threshold(threshold = method_info[["threshold"]], ms_data = ms_data, is_rt = FALSE)
            }
        }
    }
}











validate_diagnostic_frags = function(diagnostic_frags_list, diagnostic_frags_exp) {
  validate_positive_number_vector(
    x = diagnostic_frags_list, 
    name = "diagnostic_frags_list", 
    is_integers = FALSE, 
    is_named_vector = TRUE)

  if (
    !is.character(diagnostic_frags_exp) ||
    length(diagnostic_frags_exp) != 1L ||
    is.na(diagnostic_frags_exp) ||
    !nzchar(trimws(diagnostic_frags_exp))
  ) {
    stop(
      "diagnostic_frags_exp must be a single non-empty string.",
      call. = FALSE
    )
  }
  
  parsed_exp = tryCatch(
    parse(text = diagnostic_frags_exp),
    error = function(e) {
      stop(
        "Invalid diagnostic_frags_exp syntax: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (length(parsed_exp) != 1L) {
    stop(
      "diagnostic_frags_exp must contain exactly one expression.",
      call. = FALSE
    )
  }

  logical_exp_names = all.vars(parsed_exp)

  unknown_names = setdiff(logical_exp_names, names(diagnostic_frags_list))

  if (length(unknown_names) > 0L) {
    stop(
      "Unknown diagnostic fragment names in diagnostic_frags_exp: ",
      paste(unknown_names, collapse = ", "),
      call. = FALSE
    )
  }    
}



validate_mass_error_tolerance = function(mass_error_tolerance_type, mass_error_tolerance_val) {
    
    built_in_mass_error_tolerance_types = get_mass_error_tolerance_types()

    if (!(mass_error_tolerance_type %in% built_in_mass_error_tolerance_types)) {
      stop(
        "Invalid mass error tolerance type. Must be one of: ",
        paste(built_in_mass_error_tolerance_types, collapse = ", ")
      )
    }

    validate_single_nonnegative_number(
      x = mass_error_tolerance_val, 
      name = "mass_error_tolerance_val", 
      is_integers = FALSE
    )
}





validate_ms2_spectra_matching_info = function(
    diagnostic_ms2_result, 
    glycan_library, 
    mass_error_tolerance_type, 
    mass_error_tolerance_val, 
    max_possible_candidates_num
    ) {
        # validate the input parameters
        if (!inherits(diagnostic_ms2_result, "diagnostic_ms2_result")) {
            stop("diagnostic_ms2_result must be of class 'diagnostic_ms2_result'")
        }
        if (!inherits(glycan_library, "glycan_library")) {
            stop("glycan_library must be of class 'glycan_library'")
        }

        required_diagnostic_ms2_output_cols = get_required_diagnostic_ms2_output_cols()
        required_glycan_library_cols = get_required_glycan_library_cols(glycan_library = glycan_library)

        if (!all(required_diagnostic_ms2_output_cols %in% colnames(diagnostic_ms2_result))) {
            stop("diagnostic_ms2_result is missing required columns")
        }

        library_table <- glycan_library$glycan_library
        if (!is.data.frame(library_table)) {
            stop(
                "`glycan_library$glycan_library` must be a data frame.",
                call. = FALSE
            )
        }
        missing_cols <- setdiff(
            required_glycan_library_cols,
            colnames(library_table)
        )
        if (length(missing_cols) > 0L) {
            stop(
                "The glycan library table is missing required columns: ",
                paste(missing_cols, collapse = ", "),
                call. = FALSE
            )
        }

        validate_mass_error_tolerance(
            mass_error_tolerance_type = mass_error_tolerance_type,
            mass_error_tolerance_val = mass_error_tolerance_val
        )

        validate_single_positive_number(
            x = max_possible_candidates_num,
            name = "max_possible_candidates_num",
            is_integers = TRUE
        )
}





