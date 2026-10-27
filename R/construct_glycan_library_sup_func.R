resolve_monosaccharides = function(config) {
    built_in_monosaccharides = get_reference_monosaccharides()

    if (config$flags$has_custom_monosaccharide) {
        monosaccharides_weight_vector = c(built_in_monosaccharides, config$custom_monosaccharide_mass)

        custom_monosaccharide_range = list()
        for (mono_name in names(config$custom_monosaccharide_mass)) {
            custom_monosaccharide_range[[mono_name]] = c(
                config$custom_monosaccharide_min[[mono_name]],
                config$custom_monosaccharide_max[[mono_name]]
            )
        }

        monosaccharide_range = c(config$monosaccharide_range, custom_monosaccharide_range)

    } else {
        monosaccharides_weight_vector = built_in_monosaccharides

        monosaccharide_range = config$monosaccharide_range
    }
    
    # remove the monosaccharides that have zero max count
    glycan_type = config$glycan_type
    std_monosaccharide = get_glycan_type_std_monosaccharides(glycan_type)

    monosaccharides_weight_vector = monosaccharides_weight_vector[
        sapply(names(monosaccharides_weight_vector), function(mono_name) {
            mono_name %in% std_monosaccharide || 
            monosaccharide_range[[mono_name]][2] > 0
        })
    ]

    monosaccharide_range = monosaccharide_range[
        sapply(names(monosaccharide_range), function(mono_name) {
            mono_name %in% std_monosaccharide || 
            monosaccharide_range[[mono_name]][2] > 0
        })
    ]

    return(list(
        monosaccharides_weight_vector = monosaccharides_weight_vector,
        monosaccharide_range = monosaccharide_range
    ))
}


resolve_adducts = function(config) {

    if (!is.null(config$adduct_type)) {
        built_in_adducts = get_reference_adducts()
        built_in_adducts = dplyr::filter(built_in_adducts, adduct %in% config$adduct_type)
    } else {
        built_in_adducts = data.frame(
            adduct = character(),
            mass = numeric(),
            charge = numeric(),
            stringsAsFactors = FALSE
        )
    }

    if (config$flags$has_custom_adduct) {
        adducts_data = rbind(
            built_in_adducts,
            data.frame(
                adduct = names(config$custom_adduct_mass),
                mass = unname(config$custom_adduct_mass),
                charge = unname(config$custom_adduct_charge[names(config$custom_adduct_mass)]), 
                stringsAsFactors = FALSE
            )
        )
    } else {
        adducts_data = built_in_adducts
    }

    min_charge = config$min_charge_state
    max_charge = config$max_charge_state

    adducts_ranges = list()

    for (addu in adducts_data$adduct) {
        charge = adducts_data$charge[adducts_data$adduct == addu]
        adducts_ranges[[addu]] = c(0, ceiling(max_charge / charge))
    }

    return(list(
        adducts_data = adducts_data,
        adducts_ranges = adducts_ranges
    ))
}




resolve_modifications = function(config) {
    modification_weight_vector = config$custom_modification_mass

    modification_ranges = list()
    for (mod_name in names(config$custom_modification_mass)) {
        modification_ranges[[mod_name]] = c(
            config$custom_modification_min[[mod_name]],
            config$custom_modification_max[[mod_name]]
        )
    }

    return(list(
        modification_weight_vector = modification_weight_vector,
        modification_ranges = modification_ranges
    ))
}




make_mono_adduct_modif_combinations = function(monosaccharides, adducts, modifications, config) {
    monosaccharide_range = monosaccharides$monosaccharide_range

    adducts_ranges = adducts$adducts_ranges

    if (config$flags$has_custom_modification) {
        modification_ranges = modifications$modification_ranges
    } else {
        modification_ranges = list()
    }

    component_ranges = c(
        monosaccharide_range,
        adducts_ranges,
        modification_ranges
    )

    component_counts = lapply(
        component_ranges,
        function(count_range) {
            seq.int(count_range[1], count_range[2])
        }
    )

    mono_adduct_modif_combinations = expand.grid(component_counts,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE)

    return(mono_adduct_modif_combinations)
}





use_n_glycan_rules = function(combo_data, mono_name) {
    other_mono_name = setdiff(mono_name, c("Hex", "HexNAc", "dHex"))

    if ('Neu5Gc' %in% colnames(combo_data)) {
        combo_data = combo_data |> 
        dplyr::filter(
            dplyr::if_else(Neu5Ac > 0 | HexA > 0 | Neu5Gc > 0,
                         Hex >= 3 & HexNAc >= 3 & ((HexNAc-2) >= (Neu5Ac + HexA + Neu5Gc)),
                         TRUE)
          )  |>
        dplyr::filter(
            dplyr::if_else(dHex >=1 & !(Hex %in% 1:3),
                         (HexNAc + Hex -5+1) >= dHex,
                         TRUE)
        ) |>
        dplyr::filter(
            dplyr::if_else(HexNAc == 2 & Hex %in% 1:3,
                         dHex <= 1 & dplyr::if_all(dplyr::all_of(other_mono_name), ~ .x == 0),
                         TRUE)
        ) |>
        dplyr::filter(
            dplyr::if_else(HexNAc >= 3,
                         Hex >= 3,
                         TRUE)
        )

    } else {
        combo_data = combo_data |>
        dplyr::filter(
            dplyr::if_else(Neu5Ac > 0 | HexA > 0,
                         Hex >= 3 & HexNAc >= 3 & ((HexNAc-2) >= (Neu5Ac + HexA)),
                         TRUE)
          ) |>
        dplyr::filter(
            dplyr::if_else(dHex >=1 & !(Hex %in% 1:3),
                         (HexNAc + Hex -5+1) >= dHex,
                         TRUE)
          ) |>
        dplyr::filter(
            dplyr::if_else(HexNAc == 2 & Hex %in% 1:3,
                         dHex <= 1 & dplyr::if_all(dplyr::all_of(other_mono_name), ~ .x == 0),
                         TRUE)
          ) |>
        dplyr::filter(
            dplyr::if_else(HexNAc >= 3,
                         Hex >= 3,
                         TRUE)
        )

    }
    return(combo_data)
}

use_o_glycan_rules = function(combo_data) {
    if ('Neu5Gc' %in% colnames(combo_data)) {
        combo_data = combo_data |> 
        dplyr::filter(
            dplyr::if_else(Neu5Ac > 0 | Neu5Gc > 0,
                         (Hex + HexNAc + 1) >= (Neu5Ac + Neu5Gc) & Hex >= 1 & HexNAc >= 1,
                         TRUE)
        ) |>
        dplyr::filter(
            dplyr::if_else(Hex >= 0 & HexNAc > 0,
                         (HexNAc + 1) >= Hex,
                         TRUE)
        ) |>
        dplyr::filter(
            dplyr::if_else(dHex >= 0,
                         (Hex + HexNAc) >= dHex,
                         TRUE)
        ) |>
        dplyr::filter(
            dplyr::if_else(Neu5Ac > 0 | Neu5Gc > 0 | dHex >= 0,
                         (Hex + HexNAc) >= (Neu5Ac + Neu5Gc + dHex),
                         TRUE)
        ) |>
        dplyr::filter(
            dplyr::if_else(HexA > 0,
                         Hex >= HexA & Hex >= 1 & HexNAc >= 1,
                         TRUE)
        )

    } else {
      combo_data = combo_data |>
        dplyr::filter(
          dplyr::if_else(Neu5Ac > 0,
                         (Hex + HexNAc + 1) >= Neu5Ac & Hex >= 1 & HexNAc >= 1,
                         TRUE)
        ) |>
        dplyr::filter(
          dplyr::if_else(Hex >= 0 & HexNAc > 0,
                         (HexNAc + 1) >= Hex,
                         TRUE)
        ) |>
        dplyr::filter(
          dplyr::if_else(dHex >= 0,
                         (Hex + HexNAc) >= dHex,
                         TRUE)
        ) |>
        dplyr::filter(
          dplyr::if_else(Neu5Ac > 0 | dHex >= 0,
                         (Hex + HexNAc) >= (Neu5Ac + dHex),
                         TRUE)
        ) |>
        dplyr::filter(
          dplyr::if_else(HexA > 0,
                         Hex >= HexA & Hex >= 1 & HexNAc >= 1,
                         TRUE)
        )
    }
    return(combo_data)
}


use_gpi_rules = function(combo_data) {
    if (!"EtNP" %in% colnames(combo_data) || !"AHM" %in% colnames(combo_data)) {
        stop("EtNP and AHM must be defined in the `min_EtNP_num`, `max_EtNP_num`, `min_AHM_num`, and `max_AHM_num`")
    }

    if ('Neu5Gc' %in% colnames(combo_data)) {
      combo_data = combo_data |>
        dplyr::filter(
          dplyr::if_else(Neu5Ac > 0 | Neu5Gc > 0,
                         Hex >= 4,
                         TRUE)
          ) |>
        dplyr::filter(
          dplyr::if_else(EtNP > 0,
                         Hex >= EtNP,
                         TRUE)
          ) |>
        dplyr::filter(
          dplyr::if_else(
            Hex == 5,
            HexNAc == 1,
            TRUE)
          )

    } else {
      combo_data = dplyr::filter(combo_data, dplyr::if_else(Neu5Ac > 0, Hex >= 4, TRUE)) |>
        dplyr::filter(dplyr::if_else(EtNP > 0, Hex >= EtNP, TRUE)) |>
        dplyr::filter(dplyr::if_else(Hex == 5, HexNAc == 1, TRUE))
    }
    return(combo_data)
}




use_gsl_rules = function(combo_data) {
    combo_data = combo_data |>
      dplyr::filter(
        dplyr::if_else(dHex >= 1,
                       Hex >= 3,
                       TRUE)
      ) |>
      dplyr::filter(
        dplyr::if_else(HexNAc >= 3,
                       Hex >= 2,
                       TRUE)
      ) |>
      dplyr::filter(
        dplyr::if_else(dHex > 0,
                       dHex <= Hex + HexNAc - 1,
                       TRUE)
      ) |>
      dplyr::filter(
        dplyr::if_else(Hex >= 1,
                       Hex - HexNAc >= 1 & Hex - HexNAc <= 3,
                       TRUE)
      ) |>
      dplyr::filter(
        dplyr::if_else(Neu5Ac + Neu5Gc >= 2,
                       Hex >=2,
                       TRUE)
      )    
    return(combo_data)
}



use_hmo_rules = function(combo_data) {
    combo_data = combo_data |> 
        dplyr::filter(
            dplyr::if_else(dHex > 0 | Neu5Ac > 0,
                           dHex + Neu5Ac >= 0 & dHex + Neu5Ac <= 5,
                           TRUE)
        ) |> 
        dplyr::filter(
            dplyr::if_else(dHex >= 1, 
                            Hex >= 2,
                            TRUE)
        ) |> 
        dplyr::filter(
            dplyr::if_else(Neu5Ac >= 1,
                            Hex >= 2,
                            TRUE)
        ) |> 
        dplyr::filter(
            dplyr::if_else(HexNAc >= 1,
                            Hex >= HexNAc,
                            TRUE)
        ) |> 
        dplyr::filter(
            dplyr::if_else(HexNAc >= 0,
                            Hex - HexNAc <= 3,
                            TRUE)
        ) |> 
        dplyr::filter(
            dplyr::if_else(dHex > 0 | Neu5Ac > 0,
                            dHex + Neu5Ac <= Hex + HexNAc,
                            TRUE)
        )
    return(combo_data)
}









apply_glycan_library_rules = function(mono_adduct_modif_combinations, monosaccharides, adducts, modifications, config) {
    # get the names of monosaccharides, adducts, and modifications
    mono_name = names(monosaccharides$monosaccharide_range)
    adduct_name = names(adducts$adducts_ranges)
    adduct_charge = adducts$adducts_data$charge[match(adduct_name, adducts$adducts_data$adduct)
]
    # if (config$flags$has_custom_modification) {
    #     modif_name = names(modifications$modification_ranges)
    # } else {
    #     modif_name = character(0)
    # }

    # apply glycan rules
    glycan_type = config$glycan_type

    if (glycan_type == "N_glycan") {
        mono_adduct_modif_combinations = use_n_glycan_rules(combo_data = mono_adduct_modif_combinations, mono_name = mono_name)
    } else if (glycan_type == "O_glycan") {
        mono_adduct_modif_combinations = use_o_glycan_rules(combo_data = mono_adduct_modif_combinations)
    } else if (glycan_type == "GPI") {
        mono_adduct_modif_combinations = use_gpi_rules(combo_data = mono_adduct_modif_combinations)
    } else if (glycan_type == "GSL") {
        mono_adduct_modif_combinations = use_gsl_rules(combo_data = mono_adduct_modif_combinations)
    } else if (glycan_type == "HMO") {
        mono_adduct_modif_combinations = use_hmo_rules(combo_data = mono_adduct_modif_combinations)
    }

    # calculate the total number of monosaccharides 
    total_counts = rowSums(
        mono_adduct_modif_combinations[, mono_name, drop = FALSE]
    )

    mono_adduct_modif_combinations = dplyr::mutate(
        mono_adduct_modif_combinations,
        total_monosaccharides = .env$total_counts
    )

    # filter out combinations that do not meet the total monosaccharides constraints
    mono_adduct_modif_combinations = dplyr::filter(
        mono_adduct_modif_combinations,
        .data[["total_monosaccharides"]] >= .env$config$total_monosaccharide_range[1] & 
        .data[["total_monosaccharides"]] <= .env$config$total_monosaccharide_range[2]
    )

    total_charge = as.vector(
        as.matrix(
            mono_adduct_modif_combinations[, adduct_name, drop = FALSE]
        ) %*% adduct_charge
    )
    
    mono_adduct_modif_combinations = dplyr::mutate(
        mono_adduct_modif_combinations,
        charge_state = .env$total_charge
    )


    # filter out combinations that do not meet the adduct constraints
    mono_adduct_modif_combinations = dplyr::filter(
        mono_adduct_modif_combinations,
        .data[["charge_state"]] >= .env$config$min_charge_state &
        .data[["charge_state"]] <= .env$config$max_charge_state
    )


    # Allow at most one distinct non-H adduct type per combination
    non_H_adducts = setdiff(adduct_name, "H")

    if (length(non_H_adducts) > 1L) {
        non_H_adduct_type_count = rowSums(
            dplyr::select(
                mono_adduct_modif_combinations,
                dplyr::all_of(non_H_adducts)
            ) > 0
        )

        mono_adduct_modif_combinations = dplyr::filter(
            mono_adduct_modif_combinations,
            .env$non_H_adduct_type_count <= 1L
        )
    }

    # add the derivatization column based on the derivatization type specified in the config
    derivatization_type = config$derivatization_type

    mono_adduct_modif_combinations = dplyr::mutate(
        mono_adduct_modif_combinations,
        !!derivatization_type := rep(1, dplyr::n())
    )

    mono_adduct_modif_combinations = dplyr::select(
        mono_adduct_modif_combinations, 
        dplyr::all_of(mono_name), 
        dplyr::all_of(adduct_name), 
        !!derivatization_type, 
        total_monosaccharides, 
        charge_state, 
        dplyr::everything()        
    )

    return(mono_adduct_modif_combinations)
}




calculate_mono_adduct_modif_combinations_mass = function(mono_adduct_modif_combinations, monosaccharides, adducts, modifications, config) {
    monosaccharides_weight_vector = monosaccharides$monosaccharides_weight_vector
    adducts_data = adducts$adducts_data

    # calculate monosaccharide mass contribution
    monosaccharide_mass = as.vector(
        as.matrix(
            dplyr::select(
                mono_adduct_modif_combinations,
                dplyr::all_of(names(monosaccharides_weight_vector))
            )
        ) %*% monosaccharides_weight_vector
    )

    # calculate adduct mass contribution
    adduct_name = names(adducts$adducts_ranges)

    adduct_weight_vector = adducts_data$mass[match(adduct_name, adducts_data$adduct)]

    adduct_mass = as.vector(
        as.matrix(
            dplyr::select(
                mono_adduct_modif_combinations,
                dplyr::all_of(adduct_name)
            )
        ) %*% adduct_weight_vector
    )

    # calculate modification mass contribution
    modification_mass = rep(0, nrow(mono_adduct_modif_combinations))

    if (config$flags$has_custom_modification) {
        modification_weight_vector = modifications$modification_weight_vector

        modification_mass = as.vector(
            as.matrix(
                dplyr::select(
                    mono_adduct_modif_combinations,
                    dplyr::all_of(names(modification_weight_vector))
                )
            ) %*% modification_weight_vector
        )
    }

    # get the derivatization mass
    if (config$flags$has_custom_derivatization) {
        derivatization_weight = config$custom_derivatization_mass
    } else {
        build_in_derivatization_weight = get_reference_derivatizations()
        derivatization_weight = build_in_derivatization_weight[config$derivatization_type]
    }
    
    
    format_component_composition = function(data, component_names) {
        counts = as.matrix(
            dplyr::select(data, dplyr::all_of(component_names))
            )

        vapply(
            seq_len(nrow(counts)), 
            function(i) {
                quantities = counts[i, ]
                keep = quantities != 0
                
                paste0(
                    component_names[keep], quantities[keep], collapse = "")
            }, 
        character(1)
        )
    }
    
    glycan_composition = format_component_composition(
        data = mono_adduct_modif_combinations, 
        component_names = names(monosaccharides$monosaccharide_range)
    )

    adduct_type = format_component_composition(
        data = mono_adduct_modif_combinations,
        component_names = names(adducts$adducts_ranges)
    )


    # Add composition labels and calculate neutral monoisotopic mass and m/z value 
    mono_adduct_modif_combinations = dplyr::mutate(
        mono_adduct_modif_combinations, 

        glycan_composition = .env$glycan_composition, 
          
        adduct_type = .env$adduct_type,

        neutral_monoisotopic_mass = 
        .env$monosaccharide_mass +
        .env$modification_mass +
        unname(.env$derivatization_weight), 
        
        mz_value = (.data[["neutral_monoisotopic_mass"]] + .env$adduct_mass) / .data[["charge_state"]]
    )

    return(mono_adduct_modif_combinations)
}



