#---------------------------------------------
# glycan library construction main functions
#---------------------------------------------
#' Print Glycan Library Configuration Templates
#'
#' @description
#' Print a built-in configuration template for the specified glycan type.
#' The printed argument names and values can be used as a starting point
#' when calling `ConstructGlycanLibrary()`.
#'
#' @param glycan_type A single character string specifying the glycan type.
#'   Must be one of `"N_glycan"`, `"O_glycan"`, `"GPI"`, `"GSL"`, or
#'   `"HMO"`. Values are case-sensitive.
#'
#' @details
#' Each template includes:
#' \itemize{
#'   \item The glycan type and ionization mode.
#'   \item Minimum and maximum charge states and adduct types.
#'   \item The derivatization type.
#'   \item Minimum and maximum counts for the built-in components:
#'     Hex, HexNAc, dHex, Neu5Ac, HexA, Neu5Gc, Pentose, KDN, EtNP,
#'     and AHM.
#'   \item Minimum and maximum total monosaccharide counts.
#' }
#'
#' All templates specify `mode = "positive"`,
#' `adduct_type = c("H", "Na", "K")`, and
#' `derivatization_type = "ProA"`. Charge-state limits and component-count
#' limits depend on the selected glycan type.
#'
#' These templates are fixed examples and are not generated from the
#' current defaults of `ConstructGlycanLibrary()`. In particular, the
#' templates specify `"ProA"` derivatization, whereas the construction
#' function defaults to `derivatization_type = NULL`.
#'
#' This function only prints the template. It does not create a
#' configuration object, construct a glycan library, or apply the printed
#' settings to subsequent function calls. To use a template, copy its
#' arguments into a call to `ConstructGlycanLibrary()` and adjust them
#' for the experiment. Custom monosaccharide, adduct, modification, and
#' derivatization arguments are not included in the printed templates.
#'
#' An unsupported glycan type causes an error.
#'
#' @return Invisibly returns `NULL`. The configuration template is written
#'   to standard output using `cat()`.
#'
#' @seealso `ConstructGlycanLibrary()` for glycan library construction.
#'
#' @export
#'
#' @examples
#' # Print the N-glycan configuration template
#' PrintGlycanLibraryConfigsInfo("N_glycan")
#'
#' # Print the HMO configuration template
#' PrintGlycanLibraryConfigsInfo("O_glycan")
#'
#' # Capture the printed template as a character vector
#' template <- capture.output(PrintGlycanLibraryConfigsInfo("GPI"))
PrintGlycanLibraryConfigsInfo = function(glycan_type) {
  if (!glycan_type %in% c("N_glycan", "O_glycan", "GPI", "GSL", "HMO")) {
    stop("Invalid glycan type specified.")
  }

  if (glycan_type == "N_glycan") {
    cat(
      'glycan_type = "N_glycan", mode = "positive",

      min_charge_state = 1, max_charge_state = 3,
      adduct_type = c("H", "Na", "K"),

      derivatization_type = "ProA",

      min_Hex_num = 1, max_Hex_num = 12,
      min_HexNAc_num = 2, max_HexNAc_num = 10,
      min_dHex_num = 0, max_dHex_num = 4,
      min_Neu5Ac_num = 0, max_Neu5Ac_num = 4,
      min_HexA_num = 0, max_HexA_num = 2,
      min_Neu5Gc_num = 0, max_Neu5Gc_num = 4,
      min_Pentose_num = 0, max_Pentose_num = 0,
      min_KDN_num = 0, max_KDN_num = 0,
      min_EtNP_num = 0, max_EtNP_num = 0,
      min_AHM_num = 0, max_AHM_num = 0,
      min_total_monosaccharides_num = 3,
      max_total_monosaccharides_num = 22
      '
    )
  } else if (glycan_type == "O_glycan") {
    cat(
      'glycan_type = "O_glycan", mode = "positive", 

      min_charge_state = 1, max_charge_state = 2, 
      adduct_type = c("H", "Na", "K"), 

      derivatization_type = "ProA", 

      min_Hex_num = 0,      max_Hex_num = 5, 
      min_HexNAc_num = 0,   max_HexNAc_num = 5, 
      min_dHex_num = 0,     max_dHex_num = 3, 
      min_Neu5Ac_num = 0,   max_Neu5Ac_num = 3, 
      min_HexA_num = 0,     max_HexA_num = 2, 
      min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
      min_Pentose_num = 0,  max_Pentose_num = 2, 
      min_KDN_num = 0,      max_KDN_num = 0, 
      min_EtNP_num = 0, max_EtNP_num = 0,
      min_AHM_num = 0, max_AHM_num = 0,
      min_total_monosaccharides_num = 1, 
      max_total_monosaccharides_num = 10
      '
    )
  } else if (glycan_type == "GPI") {
    cat(
      'glycan_type = "GPI", mode = "positive",

      min_charge_state = 1, max_charge_state = 2, 
      adduct_type = c("H", "Na", "K"), 
      
      derivatization_type = "ProA", 

      min_Hex_num = 0,      max_Hex_num = 5, 
      min_HexNAc_num = 0,   max_HexNAc_num = 1, 
      min_dHex_num = 0,     max_dHex_num = 0, 
      min_Neu5Ac_num = 0,   max_Neu5Ac_num = 1, 
      min_HexA_num = 0,     max_HexA_num = 0, 
      min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
      min_Pentose_num = 0,  max_Pentose_num = 0, 
      min_KDN_num = 0,      max_KDN_num = 0, 
      min_EtNP_num = 0, max_EtNP_num = 3,
      min_AHM_num = 0, max_AHM_num = 1, 
      min_total_monosaccharides_num = 2, 
      max_total_monosaccharides_num = 11
      '
    )
  } else if (glycan_type == "GSL") {
    cat(
      'glycan_type = "GSL", mode = "positive",

      min_charge_state = 1, max_charge_state = 3, 
      adduct_type = c("H", "Na", "K"), 
      
      derivatization_type = "ProA", 

      min_Hex_num = 1,      max_Hex_num = 7, 
      min_HexNAc_num = 0,   max_HexNAc_num = 5, 
      min_dHex_num = 0,     max_dHex_num = 6, 
      min_Neu5Ac_num = 0,   max_Neu5Ac_num = 4, 
      min_HexA_num = 0,     max_HexA_num = 0, 
      min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
      min_Pentose_num = 0,  max_Pentose_num = 0, 
      min_KDN_num = 0,      max_KDN_num = 0, 
      min_EtNP_num = 0, max_EtNP_num = 0, 
      min_AHM_num = 0, max_AHM_num = 0,
      min_total_monosaccharides_num = 2, 
      max_total_monosaccharides_num = 16
      '
    )
  } else if (glycan_type == "HMO") {
    cat(
      'glycan_type = "HMO", mode = "positive", 
  
      min_charge_state = 1, max_charge_state = 3,   
      adduct_type = c("H", "Na", "K"),

      derivatization_type = "ProA",
      
      min_Hex_num = 2, max_Hex_num = 6,
      min_HexNAc_num = 0, max_HexNAc_num = 4,
      min_dHex_num = 0, max_dHex_num = 5,
      min_Neu5Ac_num = 0, max_Neu5Ac_num = 3,
      min_HexA_num = 0, max_HexA_num = 0,
      min_Neu5Gc_num = 0, max_Neu5Gc_num = 0,
      min_Pentose_num = 0, max_Pentose_num = 0,
      min_KDN_num = 0, max_KDN_num = 0, 
      min_EtNP_num = 0, max_EtNP_num = 0,
      min_AHM_num = 0, max_AHM_num = 0,
      min_total_monosaccharides_num = 2,
      max_total_monosaccharides_num = 13
      '
    )
  }
}




#' Construct a Glycan Composition Library
#'
#' @description
#' Enumerate glycan component, adduct, and optional modification combinations;
#' apply glycan-type-specific composition rules and count constraints; and
#' calculate neutral monoisotopic masses and theoretical m/z value. 
#' Return the resulting library together with the component
#' definitions used to construct it.
#'
#' @param glycan_type A single character string specifying the glycan type.
#'   Must be one of `"N_glycan"`, `"O_glycan"`, `"GPI"`, `"GSL"`, or
#'   `"HMO"`. Defaults to `"N_glycan"`. Changing this argument selects
#'   the composition filtering rules but does not automatically change
#'   the component-count limits.
#' @param mode A single character string, either `"positive"` or
#'   `"negative"`. Defaults to `"positive"`. In the current implementation,
#'   this argument is validated but does not change the adduct masses,
#'   charge signs, or mass-to-charge calculation. See Details.
#' @param min_charge_state,max_charge_state Positive integer-valued numeric
#'   scalars specifying the inclusive limits on total adduct charge.
#'   Defaults to `1` and `3`, respectively. The maximum must be at least
#'   the minimum.
#'
#' @param monosaccharides_additional_monoisotopic_weight_vector Optional
#'   named numeric vector of strictly positive monoisotopic mass
#'   contributions, in Da, for custom monosaccharide components.
#'   Names identify the components and become library column names.
#'   Supply together with both custom monosaccharide count vectors.
#'   Defaults to `NULL`.
#' @param monosaccharides_additional_min_num_vector Optional named numeric
#'   vector of non-negative integer minimum counts for custom monosaccharide
#'   components. Must have the same names and length as
#'   `monosaccharides_additional_monoisotopic_weight_vector`.
#'   Defaults to `NULL`.
#' @param monosaccharides_additional_max_num_vector Optional named numeric
#'   vector of non-negative integer maximum counts for custom monosaccharide
#'   components. Must have the same names and length as the corresponding
#'   mass vector. Each maximum must be at least its matching minimum.
#'   Defaults to `NULL`.
#'
#' @param adduct_type Character vector selecting built-in adducts from
#'   `"H"`, `"Na"`, `"K"`, `"Li"`, and `"Mg"`, without duplicates.
#'   Defaults to `c("H", "Na", "K")`. Built-in charges are `1` for H,
#'   Na, K, and Li, and `2` for Mg. Use `NULL` or `character(0)` to
#'   include only custom adducts; in that case, both custom adduct
#'   vectors must be supplied.
#' @param adduct_additional_monoisotopic_weight_vector Optional named
#'   numeric vector of strictly positive mass contributions, in Da,
#'   for custom adducts. Names identify the adducts and become library
#'   column names. Supply together with
#'   `adduct_additional_charge_state_vector`. Defaults to `NULL`.
#' @param adduct_additional_charge_state_vector Optional named numeric
#'   vector of positive integer charges per custom adduct. Must have
#'   the same names and length as
#'   `adduct_additional_monoisotopic_weight_vector`. Defaults to `NULL`.
#'
#' @param modification_monoisotopic_weight_vector Optional named numeric
#'   vector of strictly positive mass increments, in Da, per modification.
#'   Names identify the modifications and become library column names.
#'   Supply together with both modification count vectors.
#'   Negative mass shifts are not supported. Defaults to `NULL`.
#' @param modification_min_num_vector Optional named numeric vector of
#'   non-negative integer minimum modification counts. Must have the same
#'   names and length as `modification_monoisotopic_weight_vector`.
#'   Defaults to `NULL`.
#' @param modification_max_num_vector Optional named numeric vector of
#'   non-negative integer maximum modification counts. Must have the same
#'   names and length as `modification_monoisotopic_weight_vector`.
#'   Each maximum must be at least its matching minimum.
#'   Defaults to `NULL`.
#'
#' @param derivatization_type A single non-empty character string naming
#'   a built-in derivatization (`"ProA"`, `"AB"`, or `"PA"`) or a custom
#'   derivatization. A custom name requires
#'   `derivatization_additional_monoisotopic_weight`.
#'   Although the formal default is `NULL`, the current implementation
#'   requires an explicit name and rejects `NULL`.
#' @param derivatization_additional_monoisotopic_weight A single finite,
#'   non-negative numeric mass contribution, in Da, for a custom
#'   derivatization. Applied once per library entry. Must be `NULL`
#'   when using a built-in derivatization. Defaults to `NULL`.
#'
#' @param min_Hex_num,max_Hex_num Inclusive minimum and maximum Hex counts.
#'   Defaults to `1` and `12`.
#' @param min_HexNAc_num,max_HexNAc_num Inclusive minimum and maximum HexNAc
#'   counts. Defaults to `2` and `10`.
#' @param min_dHex_num,max_dHex_num Inclusive minimum and maximum dHex counts.
#'   Defaults to `0` and `4`.
#' @param min_Neu5Ac_num,max_Neu5Ac_num Inclusive minimum and maximum Neu5Ac
#'   counts. Defaults to `0` and `4`.
#' @param min_HexA_num,max_HexA_num Inclusive minimum and maximum HexA counts.
#'   Defaults to `0` and `2`.
#' @param min_Neu5Gc_num,max_Neu5Gc_num Inclusive minimum and maximum Neu5Gc
#'   counts. Defaults to `0` and `4`.
#' @param min_Pentose_num,max_Pentose_num Inclusive minimum and maximum
#'   Pentose counts. Both default to `0`.
#' @param min_KDN_num,max_KDN_num Inclusive minimum and maximum KDN counts.
#'   Both default to `0`.
#' @param min_EtNP_num,max_EtNP_num Inclusive minimum and maximum EtNP counts.
#'   Both default to `0`.
#' @param min_AHM_num,max_AHM_num Inclusive minimum and maximum AHM counts.
#'   Both default to `0`.
#' @param min_total_monosaccharides_num,max_total_monosaccharides_num
#'   Inclusive limits on the sum of all retained monosaccharide-component
#'   counts, including custom components and EtNP and AHM when present.
#'   Adducts, modifications, and derivatization are excluded from this sum.
#'   Defaults to `3` and `22`.
#'
#' @details
#' All count limits must be finite, non-negative integer-valued numeric
#' scalars, with each maximum greater than or equal to its minimum.
#' Total-count limits must overlap the range permitted by the individual
#' component limits.
#'
#' Custom mass and count vectors are matched by component name.
#' Names must be non-missing, non-empty, unique, and distinct from
#' built-in component names, other supplied component names, the
#' derivatization name, and output column names.
#'
#' The function generates the full Cartesian product of the supplied
#' component, adduct, and modification ranges before filtering.
#' Broad ranges can therefore require substantial memory. Use
#' `PrintGlycanLibraryConfigsInfo()` to inspect templates for each glycan
#' type and supply the desired limits explicitly.
#'
#' Glycan-type-specific composition rules are applied in addition to
#' the supplied count limits. The resulting entries describe compositions,
#' not resolved glycan structures or linkage isomers.
#'
#' Total charge is the sum of adduct counts multiplied by their charges.
#' Each entry may contain at most one distinct adduct type other than H;
#' H may coexist with that adduct type. Multiple copies of the same
#' adduct are permitted within the charge limits.
#'
#' Neutral monoisotopic mass is calculated as the sum of monosaccharide
#' component masses, modification mass increments, and one derivatization
#' mass contribution. No separate water or reducing-end mass correction
#' is added. The mass-to-charge ratio is calculated as:
#'
#' \deqn{m/z = (M + A) / z}
#'
#' where \eqn{M} is the calculated neutral monoisotopic mass, \eqn{A}
#' is the summed adduct mass contribution, and \eqn{z} is total adduct
#' charge. All supplied mass contributions should follow this convention.
#'
#' In the current implementation, selecting `mode = "negative"` leaves
#' this calculation unchanged. It does not implement deprotonation,
#' negative adduct charges, or negative-ion-specific mass corrections.
#'
#' @return
#' A list with class `"glycan_library"` containing:
#' \describe{
#'   \item{glycan_library}{
#'     A data frame with one row per retained component, adduct, and
#'     modification combination. Contains count columns for retained
#'     monosaccharide components, selected adducts, and any modifications;
#'     a column named after the derivatization with value `1`; and the
#'     derived columns described below.
#'   }
#'   \item{monosaccharides}{
#'     A list containing `monosaccharides_weight_vector`, a named vector
#'     of component masses, and `monosaccharide_range`, a named list of
#'     minimum and maximum count pairs. Components required by the
#'     glycan-type-specific rules are retained even when their maximum
#'     count is zero; other components with zero maximum counts are omitted.
#'   }
#'   \item{adducts}{
#'     A list containing `adducts_data`, a data frame with columns
#'     `adduct`, `mass`, and `charge`, and `adducts_ranges`, a named list
#'     of the count ranges used during enumeration.
#'   }
#'   \item{derivatization_type}{
#'     The supplied derivatization name.
#'   }
#'   \item{modifications}{
#'     `NULL` when modifications are not supplied; otherwise, a list
#'     containing the named mass vector `modification_weight_vector`
#'     and the named list of count ranges `modification_ranges`.
#'   }
#' }
#'
#' The derived columns in `glycan_library` are:
#' \describe{
#'   \item{total_monosaccharides}{
#'     Sum of the retained monosaccharide-component counts.
#'   }
#'   \item{charge_state}{
#'     Total positive adduct charge.
#'   }
#'   \item{glycan_composition}{
#'     Concatenated component names and non-zero counts, such as
#'     `"Hex3HexNAc2"`. Excludes adducts, modifications, and derivatization.
#'   }
#'   \item{adduct_type}{
#'     Concatenated adduct names and non-zero counts, such as `"H1Na1"`.
#'   }
#'   \item{neutral_monoisotopic_mass}{
#'     Calculated neutral monoisotopic mass in Da, including modification
#'     and derivatization contributions.
#'   }
#'   \item{mz_value}{
#'     Calculated mass-to-charge ratio.
#'   }
#' }
#'
#' Access the result table with `result$glycan_library`. The derivatization
#' name is stored in `result$derivatization_type`; the table contains a
#' count column named after that derivatization.
#'
#' @seealso
#' `PrintGlycanLibraryConfigsInfo()` for glycan-type-specific parameter
#' templates.
#'
#' @export
#'
#' @examples
#' # Construct a small N-glycan library with explicit derivatization
#' small_config <- list(
#'   glycan_type = "N_glycan",
#'   mode = "positive",
#'   min_charge_state = 1,
#'   max_charge_state = 2,
#'   adduct_type = c("H", "Na"),
#'   derivatization_type = "ProA",
#'   min_Hex_num = 3,
#'   max_Hex_num = 4,
#'   min_HexNAc_num = 2,
#'   max_HexNAc_num = 3,
#'   max_dHex_num = 1,
#'   max_Neu5Ac_num = 0,
#'   max_HexA_num = 0,
#'   max_Neu5Gc_num = 0,
#'   min_total_monosaccharides_num = 5,
#'   max_total_monosaccharides_num = 8
#' )
#'
#' result <- do.call(ConstructGlycanLibrary, small_config)
#'
#' # Inspect the library and resolved adduct definitions
#' head(result$glycan_library)
#' result$adducts$adducts_data
#'
#' # Add a hypothetical modification with a mass increment of 100 Da
#' modified_config <- small_config
#' modified_config$modification_monoisotopic_weight_vector <- c(ModX = 100)
#' modified_config$modification_min_num_vector <- c(ModX = 0)
#' modified_config$modification_max_num_vector <- c(ModX = 1)
#'
#' modified_result <- do.call(ConstructGlycanLibrary, modified_config)
#' head(modified_result$glycan_library)
ConstructGlycanLibrary <- function(
  glycan_type = "N_glycan", mode = "positive", 

  min_charge_state = 1, max_charge_state = 3,

  # define the customized monosaccharides for glycan library construction
  monosaccharides_additional_monoisotopic_weight_vector = NULL, # c('customized1' = 100, 'customized2' = 200)
  monosaccharides_additional_min_num_vector = NULL, # c('customized1' = 1, 'customized2' = 2)
  monosaccharides_additional_max_num_vector = NULL,

  adduct_type = c("H", "Na", "K"),

  # define the customized adducts for glycan library construction
  adduct_additional_monoisotopic_weight_vector = NULL, # c('customized1' = 100, 'customized2' = 200)
  adduct_additional_charge_state_vector = NULL, # c('customized1' = 1, 'customized2' = 2)

  # define the customized modifications for glycan library construction
  modification_monoisotopic_weight_vector = NULL, # c('customized1' = 100, 'customized2' = 200)
  modification_min_num_vector = NULL, # c('customized1' = 0, 'customized2' = 2)
  modification_max_num_vector = NULL, # c('customized1' = 1, 'customized2' = 2)

  # define the built-in/customized derivatization for glycan library construction
  derivatization_type = NULL,
  derivatization_additional_monoisotopic_weight = NULL,

  min_Hex_num = 1, max_Hex_num = 12,
  min_HexNAc_num = 2, max_HexNAc_num = 10,
  min_dHex_num = 0, max_dHex_num = 4,
  min_Neu5Ac_num = 0, max_Neu5Ac_num = 4,
  min_HexA_num = 0, max_HexA_num = 2,
  min_Neu5Gc_num = 0, max_Neu5Gc_num = 4,
  min_Pentose_num = 0, max_Pentose_num = 0,
  min_KDN_num = 0, max_KDN_num = 0, 

  min_EtNP_num = 0, max_EtNP_num = 0,
  min_AHM_num = 0, max_AHM_num = 0,

  min_total_monosaccharides_num = 3, max_total_monosaccharides_num = 22

) {
  config = make_glycan_library_config(
    glycan_type = glycan_type, mode = mode, 
    min_charge_state = min_charge_state, max_charge_state = max_charge_state, 
    monosaccharides_additional_monoisotopic_weight_vector = monosaccharides_additional_monoisotopic_weight_vector, 
    monosaccharides_additional_min_num_vector = monosaccharides_additional_min_num_vector,
    monosaccharides_additional_max_num_vector = monosaccharides_additional_max_num_vector, 
    adduct_type = adduct_type, 
    adduct_additional_monoisotopic_weight_vector = adduct_additional_monoisotopic_weight_vector,
    adduct_additional_charge_state_vector = adduct_additional_charge_state_vector, 
    modification_monoisotopic_weight_vector = modification_monoisotopic_weight_vector,
    modification_min_num_vector = modification_min_num_vector, modification_max_num_vector = modification_max_num_vector, 
    derivatization_type = derivatization_type, 
    derivatization_additional_monoisotopic_weight = derivatization_additional_monoisotopic_weight, 
    min_Hex_num = min_Hex_num, max_Hex_num = max_Hex_num,
    min_HexNAc_num = min_HexNAc_num, max_HexNAc_num = max_HexNAc_num,
    min_dHex_num = min_dHex_num, max_dHex_num = max_dHex_num,
    min_Neu5Ac_num = min_Neu5Ac_num, max_Neu5Ac_num = max_Neu5Ac_num,
    min_HexA_num = min_HexA_num, max_HexA_num = max_HexA_num,
    min_Neu5Gc_num = min_Neu5Gc_num, max_Neu5Gc_num = max_Neu5Gc_num,
    min_Pentose_num = min_Pentose_num, max_Pentose_num = max_Pentose_num,
    min_KDN_num = min_KDN_num, max_KDN_num = max_KDN_num,
    min_EtNP_num = min_EtNP_num, max_EtNP_num = max_EtNP_num,
    min_AHM_num = min_AHM_num, max_AHM_num = max_AHM_num, 
    min_total_monosaccharides_num = min_total_monosaccharides_num, max_total_monosaccharides_num = max_total_monosaccharides_num
    )

    config = validate_glycan_library_config(config)


    monosaccharides = resolve_monosaccharides(config)
    adducts = resolve_adducts(config)

    if (config$flags$has_custom_modification) {
        modifications = resolve_modifications(config)
    } else {
        modifications = NULL
    }

    mono_adduct_modif_combinations = make_mono_adduct_modif_combinations(
        monosaccharides = monosaccharides,
        adducts = adducts,
        modifications = modifications,
        config = config
    )

    mono_adduct_modif_combinations = apply_glycan_library_rules(
        mono_adduct_modif_combinations = mono_adduct_modif_combinations,
        monosaccharides = monosaccharides,
        adducts = adducts,
        modifications = modifications,
        config = config
    )

    mono_adduct_modif_combinations = calculate_mono_adduct_modif_combinations_mass(
        mono_adduct_modif_combinations = mono_adduct_modif_combinations,
        monosaccharides = monosaccharides,
        adducts = adducts,
        modifications = modifications,
        config = config
    )

    return(
      structure(
        list(
          glycan_library = mono_adduct_modif_combinations, 
          monosaccharides = monosaccharides,
          adducts = adducts,
          derivatization_type = derivatization_type,
          modifications = modifications
        ), 
        class = "glycan_library"
        )
    )

}


#----------------------------------------
# MS data preprocessing main functions
#----------------------------------------
# config = MSdataPreprocessingConfig(
#   ms1 = MSLevelFilter(
#     peaksCount = FilterMeanSD(2),
#     totIonCurrent = FilterQuantile(0.05, 1),
#     rtime = FilterRange(5 * 60, 50 * 60)
#   ), 
#   ms2 = MSLevelFilter(
#     peaksCount = FilterQuantile(0.1, 1),
#     totIonCurrent = FilterQuantile(0.05, 1),
#     rtime = FilterRange(5 * 60, 50 * 60)
#   )
# )

#' Preprocess MS Data Using Configurable Spectrum Filters
#'
#' @description
#' Check MS data, apply independently configured MS1 and MS2 spectrum
#' filters, and return the retained spectra together with per-variable
#' filtering records and diagnostic plots.
#'
#' @param ms_data A `Spectra` object containing both MS1 and MS2 spectra
#'   and usable MS2 precursor-charge information. Must provide `msLevel`
#'   and `spectrumId` metadata, together with numeric metadata for each
#'   configured filtering variable: `peaksCount`, `totIonCurrent`, or
#'   `rtime`.
#' @param config An object of class `"ms_data_preprocessing_config"`
#'   created by [MSdataPreprocessingConfig()]. Contains separate `ms1`
#'   and `ms2` filtering lists created by [MSLevelFilter()].
#'   Defaults to `MSdataPreprocessingConfig()`, which specifies no
#'   filtering conditions for either MS level. See Default configuration
#'   and Configuring filters below.
#'
#' @section Default configuration:
#' Omitting `config` is equivalent to supplying:
#'
#' \preformatted{
#' MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(),
#'   ms2 = MSLevelFilter()
#' )
#' }
#'
#' Both MS-level filtering lists are empty. After the input-data checks
#' succeed, all spectra are retained in their original order.
#' No automatic peak-count, intensity, or retention-time thresholds
#' are applied.
#'
#' With this configuration, `filter_records` contains empty `ms1` and
#' `ms2` lists, and `ms_variables_plots` is an empty list.
#' Input checks still require both MS1 and MS2 spectra and check MS2
#' precursor-charge information.
#'
#' To leave only one MS level unfiltered, omit its argument when calling
#' [MSdataPreprocessingConfig()] or explicitly supply `MSLevelFilter()`.
#' For example, configuring `ms1` alone does not apply its filters to
#' MS2 spectra. Use the configuration constructor rather than `NULL`
#' to request the default behavior.
#'
#' @section Configuring filters:
#' Use [MSLevelFilter()] to specify conditions for each MS level:
#' \describe{
#'   \item{peaksCount}{
#'     Number of peaks per spectrum. Supports [FilterMeanSD()],
#'     [FilterQuantile()], and [FilterRange()].
#'   }
#'   \item{totIonCurrent}{
#'     Total ion current per spectrum. Supports [FilterMeanSD()],
#'     [FilterQuantile()], and [FilterRange()].
#'   }
#'   \item{rtime}{
#'     Retention time. Supports only [FilterRange()]. Bounds must use
#'     the same units as the input retention times and lie within the
#'     observed retention-time range of the complete input object.
#'     No unit conversion is performed.
#'   }
#' }
#'
#' Arguments set to `NULL` in [MSLevelFilter()] are omitted and impose
#' no filtering condition.
#'
#' [FilterMeanSD()] determines a lower cutoff from gaps between sorted
#' values. [FilterQuantile()] determines lower and upper bounds from
#' empirical quantiles. [FilterRange()] uses fixed lower and upper
#' bounds. See the corresponding constructor documentation for
#' parameter constraints and calculation details.
#'
#' @details
#' Input-data checks run before configuration validation and filtering.
#' The function stops if either MS1 or MS2 spectra are absent, or if
#' all MS2 precursor charges are missing or all are zero. A message is
#' emitted when the precursor-charge check passes.
#'
#' Thresholds are calculated independently for each MS level and each
#' configured variable, using finite observations from the original
#' input subset for that MS level. Filters are not applied sequentially
#' to progressively reduced datasets.
#'
#' A spectrum is retained only if it passes every configured condition
#' for its MS level. Lower and upper bounds are inclusive. Non-finite
#' values (`NA`, `NaN`, `Inf`, or `-Inf`) fail any condition applied to
#' that variable. If a filtered variable has no finite observations,
#' all spectra at that MS level fail that condition.
#'
#' Invalid configurations or failures to determine filtering bounds
#' cause an error. In particular, [FilterMeanSD()] fails when finite
#' observations are available but no qualifying jump is found.
#'
#' Only MS1 and MS2 spectra are filtered. Other MS levels are retained,
#' and the original order of all retained spectra is preserved.
#' Filtering selects whole spectra; it does not alter peak arrays or
#' perform peak-level denoising.
#'
#' Diagnostic records and plot colors describe each variable's
#' individual filtering decision. A spectrum marked `"keep"` for one
#' variable can still be excluded from the final result by another
#' configured condition.
#'
#' @return A list containing:
#' \describe{
#'   \item{processed_ms_data}{
#'     A `Spectra` object containing spectra that pass all configured
#'     conditions for their MS level, together with spectra from other
#'     MS levels. Retained spectra remain in their original input order.
#'   }
#'   \item{ms_variables_plots}{
#'     A nested list of diagnostic plots indexed by MS level and
#'     configured variable, for example
#'     `ms_variables_plots$ms1$peaksCount`.
#'     Peak-count and total-ion-current entries are patchwork objects
#'     combining plots in original spectrum order and ascending value
#'     order. Retention-time entries are ggplot objects showing
#'     observations along the retention-time axis.
#'     Black indicates `"keep"` and grey indicates `"drop"` for the
#'     individual condition. Plots are returned without being
#'     automatically printed. MS levels without configured filters
#'     have no plot entry.
#'   }
#'   \item{filter_records}{
#'     A nested list with `ms1` and `ms2` elements. Each configured
#'     variable contains a data frame with one row per original
#'     spectrum at that MS level. An unfiltered MS level has an
#'     empty list.
#'   }
#' }
#'
#' Each filtering-record data frame contains:
#' \describe{
#'   \item{spectrum_index}{
#'     One-based spectrum position in the complete input object.
#'   }
#'   \item{spectrumId}{
#'     Spectrum identifier from the input metadata.
#'   }
#'   \item{variable_name}{
#'     Name of the variable evaluated by the condition.
#'   }
#'   \item{values}{
#'     Original values of that variable.
#'   }
#'   \item{status}{
#'     `"keep"` or `"drop"` for this individual condition, rather than
#'     the combined decision across all configured variables.
#'   }
#' }
#'
#' @seealso [MSdataPreprocessingConfig()], [MSLevelFilter()],
#'   [FilterMeanSD()], [FilterQuantile()], [FilterRange()].
#'
#' @export
#'
#' @examples
#' # The default configuration contains no filtering conditions
#' default_config <- MSdataPreprocessingConfig()
#' default_config$ms1
#' default_config$ms2
#'
#' # Configure MS1 and MS2 independently
#' config <- MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(
#'     peaksCount = FilterRange(20, 1000),
#'     totIonCurrent = FilterQuantile(0.05, 1)
#'   ),
#'   ms2 = MSLevelFilter(
#'     peaksCount = FilterQuantile(0.1, 1),
#'     totIonCurrent = FilterQuantile(0.05, 1)
#'   )
#' )
#'
#' # Configure MS1 only; MS2 spectra receive no filtering conditions
#' ms1_only_config <- MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(
#'     totIonCurrent = FilterQuantile(0.05, 0.95)
#'   )
#' )
#'
#' \dontrun{
#' # ms_data must be a Spectra object satisfying the input requirements
#'
#' # Run input checks and retain all spectra using the default config
#' unfiltered_result <- MSdataPreprocessing(ms_data)
#'
#' # Apply the configured MS1 and MS2 filters
#' result <- MSdataPreprocessing(ms_data, config = config)
#' processed_ms_data <- result$processed_ms_data
#'
#' # Inspect individual filtering decisions and display a diagnostic plot
#' head(result$filter_records$ms1$peaksCount)
#' print(result$ms_variables_plots$ms1$peaksCount)
#'
#' # Add retention-time filters for input with rtime stored in seconds.
#' # The observed input range must contain both 5 and 50 minutes.
#' rt_config <- MSdataPreprocessingConfig(
#'   ms1 = MSLevelFilter(rtime = FilterRange(5 * 60, 50 * 60)),
#'   ms2 = MSLevelFilter(rtime = FilterRange(5 * 60, 50 * 60))
#' )
#' rt_result <- MSdataPreprocessing(ms_data, config = rt_config)
#' }
MSdataPreprocessing = function(ms_data, config = MSdataPreprocessingConfig()) {
  # validate ms data
  ms_data_checker(ms_data)

  # validate the config
  validate_ms_data_preprocessing_config(ms_data, config)

  # perform MS data preprocessing
  preprocessed_ms_data = do_ms_data_preprocessing(ms_data = ms_data, config = config)

  # plot
  ms_variables_plots = plot_processed_ms_data_variables(record_df = preprocessed_ms_data$filter_records)

  return(
    list(
      processed_ms_data = preprocessed_ms_data$processed_ms_data,
      ms_variables_plots = ms_variables_plots, 
      filter_records = preprocessed_ms_data$filter_records
    )
  )
}





#----------------------------------------
# MS2 spectrum denoising main functions
#----------------------------------------


#' Denoise MS2 Spectra Using Configurable Intensity Thresholds
#'
#' @description
#' Determine an intensity threshold for each MS2 spectrum and remove
#' peaks at or below that threshold. Return the processed spectra,
#' per-spectrum denoising results, and the configuration used.
#'
#' @param ms_data A non-empty `Spectra` object containing at least one
#'   MS2 spectrum. Must provide `msLevel`, `spectrumId`, `acquisitionNum`,
#'   `scanIndex`, `peaksCount`, and `totIonCurrent` metadata.
#'   Spectrum identifiers should be unique across the input object.
#'   MS2 identifiers, acquisition numbers, and scan indices must be
#'   non-missing because they are checked against extracted metadata.
#'   Each MS2 peak matrix must contain numeric `mz` and `intensity`
#'   columns. Supply finite values compatible with the selected
#'   transformation and estimation method.
#'   MS1 spectra are not required; non-MS2 spectra may be included.
#' @param config An object of class `"denoising_config"` created by
#'   [DenoisingConfigAdaptive()], [DenoisingConfigSpline()],
#'   [DenoisingConfigSegmented()], [DenoisingConfigQuantile()], or
#'   [DenoisingConfigFixed()]. Defaults to `DenoisingConfigAdaptive()`.
#'   Method-specific settings are supplied through the configuration
#'   constructor rather than as additional arguments to this function.
#'
#' @section Default configuration:
#' Omitting `config` selects [DenoisingConfigAdaptive()].
#' For each eligible MS2 spectrum, the four most intense peaks are
#' selected and their m/z values are sorted.
#'
#' Spline denoising is selected when all adjacent m/z gaps are positive
#' and differ from their median by no more than a relative tolerance
#' of `0.05`. Otherwise, segmented regression is selected.
#' This is a relative spacing test, not a ppm tolerance or a check
#' for a specified isotope mass difference.
#'
#' The default intensity transformation is `log2(intensity + 1)`.
#' With the default peak-count settings, spectra containing fewer than
#' 13 peaks are assigned status `"fail"` and a threshold of zero.
#' See [DenoisingConfigAdaptive()] for all default settings.
#'
#' @section Available methods:
#' \describe{
#'   \item{Adaptive}{
#'     [DenoisingConfigAdaptive()] selects spline or segmented
#'     regression separately for each spectrum using the m/z-spacing
#'     test described above.
#'   }
#'   \item{Spline}{
#'     [DenoisingConfigSpline()] fits smoothing splines to ranked,
#'     transformed intensity values, selects a smoothing parameter
#'     using weighted evaluation criteria, and derives a threshold
#'     from the fitted first derivative.
#'   }
#'   \item{Segmented}{
#'     [DenoisingConfigSegmented()] estimates a breakpoint in ranked,
#'     transformed intensity values and maps its rounded rank to an
#'     original intensity threshold.
#'   }
#'   \item{Quantile}{
#'     [DenoisingConfigQuantile()] locates a specified quantile among
#'     transformed, distinct intensity values and maps the selected
#'     rank to an observed original intensity.
#'   }
#'   \item{Fixed}{
#'     [DenoisingConfigFixed()] applies the same absolute intensity
#'     threshold to every MS2 spectrum without transformation or
#'     spectrum-specific threshold estimation.
#'   }
#' }
#'
#' @details
#' Configuration checks run before threshold estimation. Thresholds
#' are determined independently for each MS2 spectrum.
#'
#' Spline, segmented, and quantile estimation remove duplicate original
#' intensity values before applying the configured transformation.
#' This deduplication affects estimation only; peaks with equal
#' intensities remain separate peaks during final filtering.
#'
#' Transformations are used for threshold estimation, not to replace
#' output intensities. Thresholds are expressed on the original
#' intensity scale. For each finite threshold, a peak is retained
#' only when its original intensity is finite and strictly greater
#' than the threshold. Peaks equal to the threshold are removed.
#' The m/z and intensity values of retained peaks are preserved.
#'
#' Peak filtering is registered with `Spectra::addProcessing()` on
#' the returned `Spectra` object. Spectra are not removed, even when
#' all their peaks are filtered out. The input spectrum order is
#' preserved, and non-MS2 peak arrays are left unchanged.
#'
#' For MS2 spectra with denoising status `"success"`, `peaksCount`
#' and `totIonCurrent` are recalculated from the filtered peak arrays.
#' Both become zero when no peaks remain. Other metadata, including
#' base-peak metadata, are not recalculated by this function.
#'
#' A text progress bar is displayed during threshold estimation.
#' Diagnostic plots are not created automatically; use
#' `PlotMS2DenoisingInfo()` with the original input and returned result.
#'
#' @section Insufficient peaks and estimation errors:
#' The following checks use the original number of peaks, before
#' duplicate intensity values are removed:
#' \itemize{
#'   \item Adaptive: at least
#'     `top_peaks_for_pattern + top_peaks_to_remove + 4` peaks.
#'   \item Spline: at least `top_peaks_to_remove + 4` peaks.
#'   \item Segmented: at least 10 peaks.
#' }
#'
#' Spectra failing these checks receive `denoising_status = "fail"`,
#' `denoising_method = NA`, and `ms2_noise_threshold = 0`.
#' They remain in the returned object.
#'
#' The zero threshold is still applied to failed spectra, retaining
#' only positive finite intensities. Their `peaksCount` and
#' `totIonCurrent` metadata are not recalculated, so those fields
#' may no longer match their filtered peak arrays.
#'
#' Passing a peak-count check does not guarantee enough distinct
#' intensity values for fitting. Estimation errors are not generally
#' converted into `"fail"` records and can stop the entire call.
#' Adaptive denoising does not automatically retry the alternative
#' method when the selected method raises an error.
#'
#' @return A list of class `"denoised_ms_data_ms2"` containing:
#' \describe{
#'   \item{denoising_result}{
#'     A data frame with one row per input MS2 spectrum, in the
#'     original MS2 spectrum order, and the columns described below.
#'   }
#'   \item{denoised_ms_data}{
#'     A `Spectra` object containing all input spectra in their original
#'     order, with the MS2 peak-filtering operation registered.
#'     Includes both successfully processed and failed MS2 spectra,
#'     as well as non-MS2 spectra. `peaksCount` and `totIonCurrent`
#'     are updated only for spectra marked `"success"`.
#'   }
#'   \item{config}{
#'     The supplied denoising configuration.
#'   }
#' }
#'
#' The columns of `denoising_result` are:
#' \describe{
#'   \item{ms2_spectrum_id}{
#'     Character identifier of the input MS2 spectrum.
#'   }
#'   \item{ms2_noise_threshold}{
#'     Numeric threshold on the original intensity scale.
#'     Zero is used for spectra failing the minimum peak-count check.
#'   }
#'   \item{denoising_method}{
#'     Method used: `"spline"`, `"segmented"`, `"quantile"`, or `"fixed"`.
#'     Adaptive processing records the selected method rather than
#'     `"adaptive"`. `NA` indicates a failed minimum peak-count check.
#'   }
#'   \item{denoising_status}{
#'     `"success"` when the selected threshold routine completes,
#'     or `"fail"` for an insufficient-peak result.
#'     Success does not imply that any peaks remain after filtering.
#'   }
#' }
#'
#' @seealso [DenoisingConfigAdaptive()], [DenoisingConfigSpline()],
#'   [DenoisingConfigSegmented()], [DenoisingConfigQuantile()],
#'   [DenoisingConfigFixed()], [MSdataPreprocessing()],
#'   `PlotMS2DenoisingInfo()`, `FindMS2SpectraByDiagnosticFragments()`.
#'
#' @export
#'
#' @examples
#' # Create configurations without processing spectra
#' adaptive_config <- DenoisingConfigAdaptive()
#' quantile_config <- DenoisingConfigQuantile(
#'   transform = "none",
#'   quantile_prob = 0.1
#' )
#' fixed_config <- DenoisingConfigFixed(fixed_value = 30)
#'
#' \dontrun{
#' # ms_data must be a Spectra object satisfying the input requirements
#'
#' # Use the default adaptive configuration
#' result <- MS2SpectraDenoising(ms_data)
#'
#' # Inspect thresholds, selected methods, and processing status
#' head(result$denoising_result)
#' table(result$denoising_result$denoising_status)
#'
#' # Access the returned spectra
#' denoised_spectra <- result$denoised_ms_data
#'
#' # Explicitly select a fixed threshold
#' fixed_result <- MS2SpectraDenoising(
#'   ms_data,
#'   config = fixed_config
#' )
#'
#' # Extract successfully processed MS2 spectra
#' success_ids <- with(
#'   result$denoising_result,
#'   ms2_spectrum_id[denoising_status == "success"]
#' )
#' ms2 <- Spectra::filterMsLevel(result$denoised_ms_data, 2L)
#' successful_ms2 <- ms2[as.character(ms2$spectrumId) %in% success_ids]
#'
#' # Plot a successfully processed spectrum against its original data
#' if (length(success_ids) > 0L) {
#'   plots <- PlotMS2DenoisingInfo(
#'     ms2_spectrum_id = success_ids[1],
#'     ms_data = ms_data,
#'     denoised_data = result
#'   )
#'   print(plots$ms2_spectra_plot)
#' }
#' }
MS2SpectraDenoising = function(ms_data, config = DenoisingConfigAdaptive()) {
  # validate the config
  validate_denoising_config(config)

  # determine the MS2 signal transform function
  ms2_signal_transform_func = get_transform_function(config)

  # prepare MS spectra
  ms_data_prepared = prepare_ms_spectra_data_4_ms2_denoising(ms_data, config)

  # perform MS2 spectra denoising
  denoising_results = do_ms2_spectra_denoising(
    ms_data_prepared = ms_data_prepared, 
    ms2_signal_transform_func = ms2_signal_transform_func, 
    config = config)

  # output
  denoised_data = structure(
    list(
      denoising_result = denoising_results$denoising_result,
      denoised_ms_data = denoising_results$denoised_ms_data, 
      config = config
      ), 
      class = "denoised_ms_data_ms2"
  )
  # after denoising, the denoised_ms_data contains the successfully denoised MS2 spectra and the failed denoising ones
  return(denoised_data)
}





#' Plot MS2 Spectrum Denoising Results and Diagnostics
#'
#' @description
#' Create a spectrum plot highlighting peaks at or below a saved
#' denoising threshold, together with method-specific diagnostic
#' plots for spline or segmented-regression denoising.
#'
#' @param ms2_spectrum_id A single non-missing character string identifying
#'   the MS2 spectrum to plot. Must match a `spectrumId` in `ms_data`
#'   and an `ms2_spectrum_id` in `denoised_data$denoising_result`.
#'   The corresponding `denoising_status` must be `"success"`.
#' @param ms_data The original `Spectra` object supplied to
#'   [MS2SpectraDenoising()], before applying its denoising operation.
#'   Must provide MS2 spectrum identifiers and peak matrices containing
#'   `mz` and `intensity` columns. Its complete MS2 identifier sequence,
#'   converted to character, must be identical to
#'   `denoised_data$denoising_result$ms2_spectrum_id`, including order.
#'   Supply the original object rather than a subset containing only
#'   the requested spectrum.
#' @param denoised_data An object of class `"denoised_ms_data_ms2"`
#'   returned by [MS2SpectraDenoising()] for `ms_data`. Must contain
#'   `denoising_result`, `denoised_ms_data`, and `config`.
#'   Keep the saved configuration and result table consistent with
#'   the original denoising run.
#'
#' @details
#' The spectrum plot uses the original peak data from `ms_data`.
#' For finite intensities, peaks at or below the saved
#' `ms2_noise_threshold` are colored gold (`"#E0BA48"`), and peaks
#' above it are colored black. The axes show m/z and original intensity.
#' This displays the threshold classification on the original spectrum,
#' rather than plotting separate before-and-after spectra.
#'
#' Method-specific diagnostics are selected using the recorded
#' `denoising_method` for the requested spectrum. For adaptive
#' denoising, this is the selected `"spline"` or `"segmented"` method;
#' the adaptive spacing test is not repeated.
#'
#' Spline and segmented diagnostics are recomputed from the original
#' peak data using `denoised_data$config`. They are not retrieved
#' from stored fitted models. Duplicate original intensity values
#' are removed before transformation and fitting, following the
#' corresponding denoising implementation.
#'
#' Diagnostic fits use rank index on the horizontal axis and transformed
#' intensity on the vertical axis. These indices refer to sorted
#' transformed values, not m/z positions or acquisition order.
#'
#' The spectrum plot always uses the saved threshold, whereas diagnostic
#' threshold markers use the recomputed fits. Refitting does not update
#' the saved threshold, configuration, or denoised spectra. Changed data,
#' settings, or fitting behavior can therefore produce diagnostic fits
#' that differ from those used in the original run.
#'
#' Quantile and fixed-threshold methods produce the spectrum plot
#' without additional fitting diagnostics.
#'
#' The function returns plot objects without automatically printing
#' or saving them. Use `print()` to display individual plots.
#'
#' @section Input checks and fitting errors:
#' The requested spectrum must have denoising status `"success"`.
#' Spectra marked `"fail"` are rejected.
#'
#' Validation compares the complete MS2 identifier sequence in
#' `ms_data` with the result table, not just the requested identifier.
#' Reordering or subsetting the MS2 input without a corresponding
#' denoising run causes validation to fail. Identifier matching does
#' not verify that peak values are unchanged.
#'
#' Errors during diagnostic refitting can stop the call even when
#' the saved denoising status is `"success"`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{ms2_spectra_plot}{
#'     A ggplot object showing the original MS2 spectrum, with peak
#'     colors determined by the saved intensity threshold.
#'   }
#'   \item{denoising_info_plot}{
#'     A named list of method-specific plot objects for spline or
#'     segmented denoising, or `NULL` for quantile and fixed methods.
#'     Its contents are described below.
#'   }
#' }
#'
#' For spline denoising, `denoising_info_plot` contains:
#' \describe{
#'   \item{spar_metric_plots}{
#'     A patchwork object combining four panels in one row:
#'     RMSE, CV, D1, and D2 versus the candidate smoothing parameter
#'     `spar`. The panels show the unnormalized metric values;
#'     non-finite values are omitted. A dashed vertical line marks
#'     the selected smoothing parameter.
#'   }
#'   \item{optimal_spline_plot}{
#'     A ggplot object showing sorted transformed intensity values
#'     and the selected smoothing spline. A dashed vertical line
#'     marks the threshold rank, with different line colors on
#'     either side.
#'   }
#' }
#'
#' For segmented denoising, `denoising_info_plot` contains:
#' \describe{
#'   \item{segmented_regression_plot}{
#'     A ggplot object showing sorted transformed intensity values,
#'     the initial linear fit, and the auxiliary fits below and above
#'     the estimated breakpoint. A dashed vertical line marks the
#'     rounded breakpoint rank used to determine the threshold.
#'   }
#' }
#'
#' @seealso [MS2SpectraDenoising()], [DenoisingConfigAdaptive()],
#'   [DenoisingConfigSpline()], [DenoisingConfigSegmented()],
#'   [DenoisingConfigQuantile()], [DenoisingConfigFixed()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # ms_data is the original Spectra object used for denoising
#' denoised_data <- MS2SpectraDenoising(
#'   ms_data,
#'   config = DenoisingConfigAdaptive()
#' )
#'
#' # Select a successfully processed MS2 spectrum
#' success_ids <- with(
#'   denoised_data$denoising_result,
#'   ms2_spectrum_id[denoising_status == "success"]
#' )
#'
#' if (length(success_ids) > 0L) {
#'   plots <- PlotMS2DenoisingInfo(
#'     ms2_spectrum_id = success_ids[1],
#'     ms_data = ms_data,
#'     denoised_data = denoised_data
#'   )
#'
#'   # Display the original spectrum with threshold-based peak colors
#'   print(plots$ms2_spectra_plot)
#'
#'   # Display each available method-specific diagnostic
#'   if (!is.null(plots$denoising_info_plot)) {
#'     for (diagnostic_plot in plots$denoising_info_plot) {
#'       print(diagnostic_plot)
#'     }
#'   }
#' }
#' }
PlotMS2DenoisingInfo = function(ms2_spectrum_id, ms_data, denoised_data) {
  # check the input 
  validate_plot_ms2_denoising_info(
    ms2_spectrum_id = ms2_spectrum_id, 
    ms_data = ms_data, 
    denoised_data = denoised_data
  )

  # determine the MS2 signal transform function
  ms2_signal_transform_func = get_transform_function(denoised_data$config)

  # plot ms2 spectra
  ms2_spectra_plot = plot_ms2_spectra(
    ms2_spectrum_id = ms2_spectrum_id,
    ms_data = ms_data,
    denoised_data = denoised_data
  )

  # plot denoising info for the adaptive, spline, segmented
  denoising_info_plot = plot_ms2_denoising_info(
    ms2_spectrum_id = ms2_spectrum_id,
    ms_data = ms_data,
    denoised_data = denoised_data, 
    ms2_signal_transform_func = ms2_signal_transform_func
  )

  return(
    list(
      ms2_spectra_plot = ms2_spectra_plot,
      denoising_info_plot = denoising_info_plot
    )
  )
}



#-----------------------------------------------------------
# find MS2 spectrum by diagnostic fragments main functions
#-----------------------------------------------------------
# logical expression-based filtering 

#' Find MS2 Spectra Matching Diagnostic Fragment Conditions
#'
#' @description
#' Screen successfully denoised MS2 spectra using a logical expression
#' of diagnostic fragment presence. For each matching spectrum, return
#' its metadata together with the identifier and retention time of the
#' nearest preceding MS1 spectrum.
#'
#' @param denoised_data An object of class `"denoised_ms_data_ms2"`
#'   returned by [MS2SpectraDenoising()]. Must contain `denoising_result`
#'   and `denoised_ms_data`. Only spectra marked with
#'   `denoising_status = "success"` are screened. Supply an object
#'   containing at least one successfully processed MS2 spectrum.
#' @param ms_data_raw A `Spectra` object containing the corresponding
#'   original MS data, including MS1 spectra. Every successfully
#'   denoised MS2 identifier must occur in its `spectrumId` metadata.
#'   MS1 spectra must provide identifiers and valid retention times.
#'   Each matching MS2 spectrum must have at least one MS1 spectrum
#'   with a strictly earlier retention time.
#' @param diagnostic_frags_list A non-empty named numeric vector of
#'   finite positive target fragment m/z values. Names identify the
#'   fragments used in `diagnostic_frags_exp` and must be unique and
#'   non-empty. Defaults to [DiagnosticFragmentsConfigDefault()].
#'   Despite the argument name, supply a numeric vector rather than
#'   an R list.
#' @param diagnostic_frags_exp A single non-empty character string
#'   containing one R expression that evaluates to a single non-missing
#'   logical value for each spectrum. Referenced fragment names must
#'   occur in `diagnostic_frags_list` and are case-sensitive.
#'   Combine presence conditions using operators such as `&`, `|`,
#'   and `!`, with parentheses as needed. For example,
#'   `"HexNAc & (HexNAc_ProA | Hex)"`. This argument has no default.
#' @param mass_error_tolerance_type A single character string specifying
#'   the fragment matching tolerance: `"ppm_tolerance"` for relative
#'   error in parts per million, or `"absolute_tolerance"` for an
#'   absolute m/z difference. This argument has no default.
#' @param mass_error_tolerance_val A single finite, non-negative numeric
#'   tolerance. Express in ppm for `"ppm_tolerance"`, or in the same
#'   m/z units as the peak and target values for `"absolute_tolerance"`.
#'   Matching includes the tolerance boundary. This argument has
#'   no default.
#'
#' @section Diagnostic fragment matching:
#' Fragment presence is evaluated using the peak m/z values in
#' `denoised_data$denoised_ms_data`, after denoising. Failed spectra
#' are excluded. The original MS2 peak arrays in `ms_data_raw` are
#' not used for fragment detection.
#'
#' For an observed peak m/z value `observed` and a target fragment
#' m/z value `target`, a match satisfies:
#' \itemize{
#'   \item Relative tolerance:
#'     `abs((observed - target) / target) * 1e6 <= mass_error_tolerance_val`.
#'   \item Absolute tolerance:
#'     `abs(observed - target) <= mass_error_tolerance_val`.
#' }
#'
#' A fragment is present if at least one peak matches its target.
#' Only fragments referenced in `diagnostic_frags_exp` are evaluated.
#' Including additional fragments in `diagnostic_frags_list` does
#' not require them to be present.
#'
#' Each fragment name in the expression represents a logical presence
#' value. For example, `"HexNAc & HexNAc_ProA"` requires both fragments,
#' `"HexNAc | Hex"` requires either fragment, and
#' `"HexNAc & !Bisecting"` requires HexNAc without a detected
#' Bisecting fragment.
#'
#' No additional intensity cutoff or fragment charge correction is
#' applied. Target values are used directly as m/z values.
#' Fragment tests are independent: the same observed peak can satisfy
#' multiple targets when their tolerance windows overlap.
#' Supply finite peak m/z values; missing values are not explicitly
#' removed before presence evaluation.
#'
#' @section Association with MS1 spectra:
#' For each matching MS2 spectrum, the function selects the MS1
#' spectrum in `ms_data_raw` with the smallest positive retention-time
#' difference from that MS2 spectrum.
#'
#' MS1 spectra at the same retention time as the MS2 spectrum are
#' excluded. If several MS1 spectra share the closest earlier time,
#' the first in the input MS1 sequence is selected. No maximum
#' retention-time gap is imposed.
#'
#' Association uses retention time only; it does not inspect an
#' explicit precursor-scan link. Retention times in both objects
#' must use the same units. If a matching MS2 spectrum has no earlier
#' MS1 spectrum, the function stops with an error.
#'
#' @details
#' Fragment definitions, expression syntax, referenced fragment names,
#' and mass-error tolerance are checked before screening. The expression
#' is parsed and evaluated as R code against fragment presence values.
#'
#' Successfully denoised MS2 identifiers must exist in both
#' `ms_data_raw` and `denoised_data$denoised_ms_data`. Spectrum identifiers
#' should uniquely identify scans across the supplied objects.
#' MS2 retention times are taken from the denoised object.
#'
#' A text progress bar is displayed during screening. Matching spectra
#' are returned in their order in the denoised input. The result
#' contains metadata rather than peak arrays, fragment-presence flags,
#' matched peak m/z values, or fragment mass errors.
#'
#' @section No matching spectra:
#' If screening completes without a match, the function emits
#' `"No MS2 spectrum was identified"` and returns `NULL`.
#'
#' @return
#' When matches are found, a data frame with additional class
#' `"diagnostic_ms2_result"` and one row per matching MS2 spectrum.
#' Contains the following columns:
#' \describe{
#'   \item{ms2_spectrum_id}{
#'     Identifier of the matching MS2 spectrum.
#'   }
#'   \item{ms2_precursor_charge}{
#'     MS2 precursor charge from `precursorCharge`.
#'   }
#'   \item{ms2_precursor_mz}{
#'     MS2 precursor m/z from `precursorMz`.
#'   }
#'   \item{ms2_retention_time}{
#'     MS2 retention time from `rtime`.
#'   }
#'   \item{ms1_spectrum_id}{
#'     Identifier of the nearest strictly preceding MS1 spectrum.
#'   }
#'   \item{ms1_retention_time}{
#'     Retention time of the associated MS1 spectrum.
#'   }
#'   \item{ms2_isolation_window_target_mz}{
#'     MS2 isolation-window target m/z from `isolationWindowTargetMz`.
#'   }
#'   \item{ms2_peaks_count}{
#'     Peak count stored in the denoised MS2 `peaksCount` metadata.
#'   }
#'   \item{ms2_tot_ion_current}{
#'     Total ion current stored in the denoised MS2 `totIonCurrent`
#'     metadata.
#'   }
#'   \item{ms2_base_peak_mz}{
#'     MS2 base-peak m/z from `basePeakMZ`.
#'   }
#'   \item{ms2_base_peak_intensity}{
#'     MS2 base-peak intensity from `basePeakIntensity`.
#'   }
#'   \item{ms2_collision_energy}{
#'     MS2 collision energy from `collisionEnergy`.
#'   }
#'   \item{ms2_ionisation_energy}{
#'     MS2 ionisation energy from `ionisationEnergy`.
#'   }
#'   \item{ms2_precursor_intensity}{
#'     MS2 precursor intensity from `precursorIntensity`.
#'   }
#' }
#' Returns `NULL` if screening completes without a matching spectrum.
#'
#' MS2 metadata are copied from the denoised object, and MS1 metadata
#' from `ms_data_raw`. Missing optional metadata fields are filled
#' with `NA`. This function does not recalculate metadata; in particular,
#' base-peak fields may still describe the spectrum before denoising.
#' Retention-time and other metadata units are preserved.
#'
#' The result can be supplied as `diagnostic_ms2_result` to
#' `MatchMS2SpectraToGlycanCompositions()`.
#'
#' @seealso [DiagnosticFragmentsConfigDefault()],
#'   [MS2SpectraDenoising()], [PlotMS2DenoisingInfo()],
#'   `MatchMS2SpectraToGlycanCompositions()`.
#'
#' @export
#'
#' @examples
#' # Inspect the default target m/z values
#' fragments <- DiagnosticFragmentsConfigDefault()
#' fragments[c("HexNAc", "HexNAc_ProA", "Hex")]
#'
#' # Require HexNAc and at least one additional diagnostic fragment
#' fragment_expression <- "HexNAc & (HexNAc_ProA | Hex)"
#'
#' \dontrun{
#' # ms_data_raw contains the original MS1 and MS2 spectra
#' denoised_data <- MS2SpectraDenoising(
#'   ms_data_raw,
#'   config = DenoisingConfigFixed(fixed_value = 30)
#' )
#'
#' # Use default fragment definitions and a 10 ppm matching tolerance
#' matches <- FindMS2SpectraByDiagnosticFragments(
#'   denoised_data = denoised_data,
#'   ms_data_raw = ms_data_raw,
#'   diagnostic_frags_exp = fragment_expression,
#'   mass_error_tolerance_type = "ppm_tolerance",
#'   mass_error_tolerance_val = 10
#' )
#'
#' head(matches)
#'
#' # Use a subset of targets and an absolute m/z tolerance
#' selected_fragments <- fragments[c("HexNAc", "HexNAc_ProA")]
#' matches_absolute <- FindMS2SpectraByDiagnosticFragments(
#'   denoised_data = denoised_data,
#'   ms_data_raw = ms_data_raw,
#'   diagnostic_frags_list = selected_fragments,
#'   diagnostic_frags_exp = "HexNAc & HexNAc_ProA",
#'   mass_error_tolerance_type = "absolute_tolerance",
#'   mass_error_tolerance_val = 0.01
#' )
#' }
FindMS2SpectraByDiagnosticFragments = function(
  denoised_data, 
  ms_data_raw, 
  diagnostic_frags_list = DiagnosticFragmentsConfigDefault(), 
  diagnostic_frags_exp, 
  mass_error_tolerance_type, 
  mass_error_tolerance_val) {
    # validate the input parameters 
    validate_diagnostic_frags(
      diagnostic_frags_list = diagnostic_frags_list, 
      diagnostic_frags_exp = diagnostic_frags_exp
    )

    validate_mass_error_tolerance(
      mass_error_tolerance_type = mass_error_tolerance_type,
      mass_error_tolerance_val = mass_error_tolerance_val
    )

    # prepare the successfully denoised MS2 spectra from the denoised data    
    denoised_ms_data_ms2 = prepare_denoised_ms_data_ms2(
      denoised_data = denoised_data, 
      ms_data_raw = ms_data_raw
    )

    target_spectra = find_ms2_spectra_by_diagnostic_fragments(
      denoised_ms_data_ms2 = denoised_ms_data_ms2,
      ms_data_raw = ms_data_raw,
      diagnostic_frags_list = diagnostic_frags_list,
      diagnostic_frags_exp = diagnostic_frags_exp,
      mass_error_tolerance_type = mass_error_tolerance_type,
      mass_error_tolerance_val = mass_error_tolerance_val
    )
    # Return early when no matching MS2 spectra are found
    if (is.null(target_spectra)) {
      return(NULL)
    }

    return(
      structure(
        target_spectra, 
        class = unique(c("diagnostic_ms2_result", class(target_spectra)))
      )
    )

}



#' Match MS2 Precursors to Candidate Glycan Compositions
#'
#' @description
#' Match diagnostic MS2 precursor m/z values and charge states against
#' a constructed glycan library. Rank matching library entries by
#' absolute mass error, retain a configurable number of candidates
#' per spectrum, and attach the corresponding MS2 and MS1 metadata.
#'
#' @param diagnostic_ms2_result A data frame inheriting from
#'   `"diagnostic_ms2_result"`, typically returned by
#'   [FindMS2SpectraByDiagnosticFragments()]. Must contain
#'   `ms2_spectrum_id`, `ms2_precursor_charge`, `ms2_precursor_mz`,
#'   `ms2_retention_time`, `ms1_spectrum_id`, and `ms1_retention_time`.
#'   Supply one row per MS2 spectrum with unique, non-missing spectrum
#'   identifiers and valid numeric precursor m/z and charge values.
#' @param glycan_library An object of class `"glycan_library"` returned
#'   by [ConstructGlycanLibrary()]. Supply the complete object.
#'   Its `glycan_library` element must be a data frame containing
#'   `glycan_composition`, `adduct_type`, `neutral_monoisotopic_mass`,
#'   `mz_value`, `charge_state`, and the count column named by the
#'   object's `derivatization_type` element, such as `"ProA"`.
#'   Each table row represents one candidate library entry.
#' @param mass_error_tolerance_type A single character string specifying
#'   the matching tolerance: `"ppm_tolerance"` for relative error in
#'   parts per million, or `"absolute_tolerance"` for an absolute
#'   m/z difference.
#' @param mass_error_tolerance_val A single finite, non-negative numeric
#'   tolerance. Express in ppm for `"ppm_tolerance"`, or in the same
#'   m/z units as the observed and theoretical values for
#'   `"absolute_tolerance"`. Candidates exactly on the tolerance
#'   boundary are included.
#' @param max_possible_candidates_num A single finite, positive
#'   integer-valued numeric value specifying the maximum number of
#'   candidate library rows retained per input MS2 row.
#'   Defaults to `3`. Candidates are selected after sorting by
#'   increasing mass error.
#'
#' @section Glycan library input:
#' The function validates the complete library object and performs
#' matching against its internal `glycan_library` data frame.
#' The outer `derivatization_type` metadata identifies the required
#' derivatization count column in that table.
#'
#' For example, an object constructed with
#' `derivatization_type = "ProA"` must contain a `ProA` column in its
#' internal table. A separate table column named `derivatization_type`
#' is not required or added by this function.
#'
#' The component definitions stored in the outer object remain
#' available for subsequent analysis, but are not included in the
#' returned matching table.
#'
#' @section Matching and ranking:
#' For each input MS2 spectrum, library entries are first restricted
#' to those whose `charge_state` exactly equals `ms2_precursor_charge`.
#' Mass error is then calculated between the observed precursor m/z
#' and each remaining theoretical `mz_value`.
#'
#' For an observed precursor m/z value `observed` and a theoretical
#' library m/z value `theoretical`, the error is:
#' \itemize{
#'   \item Relative error:
#'     `abs(observed - theoretical) / theoretical * 1e6`.
#'   \item Absolute error:
#'     `abs(observed - theoretical)`.
#' }
#'
#' Entries with error less than or equal to `mass_error_tolerance_val`
#' are sorted by increasing error. The first
#' `max_possible_candidates_num` rows are retained, or all matching
#' rows if fewer are available.
#'
#' Equal-error candidates retain their library order. Ties at the
#' candidate limit are not expanded. The limit applies to library
#' rows rather than distinct `glycan_composition` labels. Different
#' adduct or modification combinations and duplicate entries are
#' not merged.
#'
#' Matching and ranking use precursor charge and m/z only.
#' Fragment identities, retention times, isotope patterns, and
#' fragment intensities do not contribute to the ranking.
#'
#' @details
#' Input classes, required columns, tolerance settings, and the
#' candidate limit are checked before matching. The internal library
#' table must be a data frame.
#'
#' Candidate rows are collected in input MS2 row order, with candidates
#' within each spectrum ordered by increasing mass error. The six
#' required diagnostic-result columns are then joined by
#' `ms2_spectrum_id`.
#'
#' Spectrum identifiers in `diagnostic_ms2_result` should be unique.
#' Uniqueness is not checked by the current validator; duplicate
#' identifiers can multiply output rows during the metadata join.
#'
#' If an individual spectrum has no candidates, a message identifies
#' that spectrum and it contributes no output rows. If no spectrum
#' has any candidates, the function stops with
#' `"No candidates found for any MS2 spectrum."`.
#' A zero-row diagnostic input also reaches this no-candidate error.
#'
#' @return
#' A data frame containing one row per retained spectrum-library
#' match. Columns from `glycan_library$glycan_library` are retained,
#' including component counts, the derivatization count column,
#' composition and adduct labels, theoretical masses, and charge states.
#' The following columns are added:
#' \describe{
#'   \item{mass_error}{
#'     Non-negative precursor m/z error, expressed in ppm or absolute
#'     m/z units according to `mass_error_tolerance_type`.
#'     This is an unsigned error.
#'   }
#'   \item{ms2_spectrum_id}{
#'     Identifier of the matched MS2 spectrum.
#'   }
#'   \item{ms2_precursor_charge}{
#'     Observed precursor charge from `diagnostic_ms2_result`.
#'   }
#'   \item{ms2_precursor_mz}{
#'     Observed precursor m/z from `diagnostic_ms2_result`.
#'   }
#'   \item{ms2_retention_time}{
#'     MS2 retention time from `diagnostic_ms2_result`.
#'   }
#'   \item{ms1_spectrum_id}{
#'     Associated MS1 spectrum identifier from `diagnostic_ms2_result`.
#'   }
#'   \item{ms1_retention_time}{
#'     Associated MS1 retention time from `diagnostic_ms2_result`.
#'   }
#' }
#'
#' Additional diagnostic-result metadata columns are not included.
#' The tolerance type is not stored as a separate output column,
#' and no explicit candidate-rank column is added.
#' No dedicated matching-result class is assigned.
#'
#' The result can be supplied as `matching_results` to
#' `ValidateGlycanCompositionByIsotopePattern()`, together with
#' the original constructed library object.
#'
#' @seealso [FindMS2SpectraByDiagnosticFragments()],
#'   [ConstructGlycanLibrary()],
#'   `ValidateGlycanCompositionByIsotopePattern()`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Construct a small N-glycan library
#' constructed <- ConstructGlycanLibrary(
#'   glycan_type = "N_glycan",
#'   derivatization_type = "ProA",
#'   min_Hex_num = 3,
#'   max_Hex_num = 4,
#'   min_HexNAc_num = 2,
#'   max_HexNAc_num = 3,
#'   max_dHex_num = 1,
#'   max_Neu5Ac_num = 0,
#'   max_HexA_num = 0,
#'   max_Neu5Gc_num = 0,
#'   min_total_monosaccharides_num = 5,
#'   max_total_monosaccharides_num = 8
#' )
#'
#' # diagnostic_result is a non-empty result returned by
#' # FindMS2SpectraByDiagnosticFragments().
#' # Pass the complete constructed library object.
#' candidates <- MatchMS2SpectraToGlycanCompositions(
#'   diagnostic_ms2_result = diagnostic_result,
#'   glycan_library = constructed,
#'   mass_error_tolerance_type = "ppm_tolerance",
#'   mass_error_tolerance_val = 10,
#'   max_possible_candidates_num = 3
#' )
#'
#' head(candidates)
#'
#' # Retain only the closest candidate within an absolute m/z tolerance
#' closest_candidates <- MatchMS2SpectraToGlycanCompositions(
#'   diagnostic_ms2_result = diagnostic_result,
#'   glycan_library = constructed,
#'   mass_error_tolerance_type = "absolute_tolerance",
#'   mass_error_tolerance_val = 0.01,
#'   max_possible_candidates_num = 1
#' )
#' }
MatchMS2SpectraToGlycanCompositions = function(
  diagnostic_ms2_result, 
  glycan_library, 
  mass_error_tolerance_type, 
  mass_error_tolerance_val, 
  max_possible_candidates_num = 3
  ) 
  {
  # validate the input parameters
  validate_ms2_spectra_matching_info(
    diagnostic_ms2_result = diagnostic_ms2_result,
    glycan_library = glycan_library,
    mass_error_tolerance_type = mass_error_tolerance_type,
    mass_error_tolerance_val = mass_error_tolerance_val,
    max_possible_candidates_num = max_possible_candidates_num
  )

  matching_results = match_ms2_spectra_to_glycan_library(
    diagnostic_ms2_result = diagnostic_ms2_result,
    glycan_library = glycan_library,
    mass_error_tolerance_type = mass_error_tolerance_type,
    mass_error_tolerance_val = mass_error_tolerance_val,
    max_possible_candidates_num = max_possible_candidates_num
  )
  return(matching_results)
}













#' Select Glycan Candidates by MS1 Isotope Pattern Similarity
#'
#' @description
#' Resolve multiple glycan candidates for an MS2 spectrum by comparing
#' their theoretical isotope patterns with the associated experimental
#' MS1 spectrum. Retain the candidate with the highest cosine similarity.
#' Spectra with only one candidate are returned without isotope scoring.
#'
#' @param matching_results A non-empty data frame of candidate matches,
#'   typically returned by [MatchMS2SpectraToGlycanCompositions()].
#'   Must contain `ms1_spectrum_id` and `ms2_spectrum_id`.
#'   For spectra with multiple candidates, also requires
#'   `ms2_precursor_mz`, `charge_state`, and count columns for every
#'   component identified by `constructed_glycan_library`.
#'   Each such MS2 group must have exactly one associated MS1 identifier
#'   and one observed precursor m/z value.
#' @param molecular_formula A named character vector of component
#'   chemical formulas to supplement the built-in definitions.
#'   Names must match component-count columns and component names in
#'   `constructed_glycan_library`. For example,
#'   `c(CustomSugar = "C6H10O5")`.
#'   Supply formulas for custom monosaccharides, adducts, modifications,
#'   and derivatizations used by the library.
#'   Pass `NULL` explicitly when the built-in definitions are sufficient;
#'   this argument has no default.
#' @param constructed_glycan_library The complete `"glycan_library"`
#'   object returned by [ConstructGlycanLibrary()] and used to generate
#'   `matching_results`. Its monosaccharide, adduct, modification, and
#'   derivatization metadata identify the component-count columns used
#'   to assemble candidate formulas.
#' @param bin_width A single finite, positive numeric spacing between
#'   m/z grid centers used to compare experimental and theoretical
#'   patterns. Defaults to `0.3`. Express in the same m/z units as the
#'   spectra. This is grid spacing, not a ppm tolerance.
#' @param threshold_iso_probalility A numeric probability cutoff passed
#'   as `threshold` to [enviPat::isopattern()] with `rel_to = 2`.
#'   Defaults to `0.01`. This is an absolute probability cutoff,
#'   not a percentage relative to the most abundant isotope peak.
#'   A value of zero requests calculation without probability pruning.
#'   The argument name retains the spelling used in the function
#'   signature.
#' @param ms_data_raw A `Spectra` object containing the original MS1
#'   spectra referenced by `matching_results`. Spectrum identifiers
#'   must uniquely identify scans. Each MS1 scan used for scoring
#'   must provide numeric `mz` and `intensity` peak columns.
#'   Supply finite, non-negative intensities and extraction windows
#'   with positive total intensity.
#' @param ms1_window_left,ms1_window_right Finite, non-negative numeric
#'   distances below and above the observed MS2 precursor m/z defining
#'   the experimental MS1 extraction window. Supply values whose sum
#'   is positive. The inclusive interval is
#'   `[ms2_precursor_mz - ms1_window_left,
#'   ms2_precursor_mz + ms1_window_right]`.
#'   Distances use m/z units and are not derived from isolation-window
#'   metadata. Both arguments must be supplied.
#'
#' @section Candidate formula construction:
#' Component names are obtained from the constructed library's
#' monosaccharide ranges, adduct ranges, optional modification ranges,
#' and derivatization name.
#'
#' User-supplied formulas are combined with built-in formulas for
#' monosaccharides, adducts, and derivatizations. User entries precede
#' the built-in definitions, so a supplied formula with the same name
#' takes precedence. Use unique names within `molecular_formula`.
#'
#' For each candidate, component formulas are multiplied by their
#' counts using [enviPat::multiform()] and combined using
#' [enviPat::mergeform()]. Components with zero counts are skipped.
#' The assembled formula is processed by [enviPat::check_chemform()].
#'
#' Formulas describe each component's contribution to the complete
#' candidate formula. No additional water or reducing-end correction
#' is added. Custom formulas should be consistent with the mass
#' contributions used to construct the glycan library; this consistency
#' is not checked by the function.
#'
#' The theoretical isotope pattern is calculated using the candidate's
#' `charge_state`, the isotope data supplied by enviPat, and
#' `rel_to = 2`. Formula-based isotope calculation does not use the
#' stored `neutral_monoisotopic_mass` to construct the pattern.
#'
#' @section Experimental window and similarity:
#' Candidates are grouped by `ms2_spectrum_id`. For each group with
#' multiple candidates, the associated MS1 scan is retrieved by
#' `ms1_spectrum_id`.
#'
#' Experimental peaks are extracted within the inclusive window around
#' the observed `ms2_precursor_mz`. Their intensities are divided by
#' the total intensity within that window.
#'
#' For each candidate, an evenly spaced grid is generated from the
#' minimum to the maximum m/z across its theoretical pattern and the
#' extracted experimental peaks, using `bin_width` as the increment.
#' Grid limits are therefore candidate-specific.
#'
#' Each experimental and theoretical peak is assigned to its nearest
#' grid center. Intensities or abundances assigned to the same center
#' are summed. No separate maximum peak-to-center distance is imposed.
#' Theoretical isotope peaks are not cropped to the experimental
#' extraction window before comparison.
#'
#' Cosine similarity is calculated between the resulting experimental
#' and theoretical vectors:
#'
#' \deqn{
#'   \mathrm{score} =
#'   \frac{\sum_i E_i T_i}
#'   {\sqrt{\sum_i E_i^2}\sqrt{\sum_i T_i^2}}
#' }
#'
#' where \eqn{E_i} and \eqn{T_i} are the experimental and theoretical
#' values at grid center \eqn{i}. For non-negative, non-zero finite
#' vectors, scores lie between zero and one, with larger values
#' indicating greater similarity.
#'
#' @section Candidate selection:
#' A group containing only one candidate is copied directly to the
#' output. Its MS1 window and theoretical isotope pattern are not
#' evaluated, and no new score is calculated.
#'
#' For a group with multiple candidates, the row with the largest
#' score is retained. Ties are resolved by keeping the first tied
#' candidate in the input order. No minimum acceptable similarity
#' threshold is applied.
#'
#' Selection operates on the candidates supplied in `matching_results`.
#' Candidates excluded during precursor-mass matching cannot be
#' recovered by this function. Retaining a single candidate does not
#' establish a glycan structure or provide a statistical confidence
#' estimate.
#'
#' @details
#' The function explicitly checks for `ms1_spectrum_id` and
#' `ms2_spectrum_id` columns. For groups with multiple candidates,
#' it also checks that the associated MS1 identifier and observed
#' precursor m/z are each unique.
#'
#' Other input requirements are not comprehensively validated.
#' Missing component formulas or count columns, unavailable MS1
#' scans, invalid grid settings, or failed isotope calculations
#' can cause errors.
#'
#' Empty extraction windows, zero total intensity, and zero-norm
#' comparison vectors are not explicitly handled and can produce
#' errors or non-finite scores. Provide suitable MS1 data and
#' window settings before calling the function.
#'
#' A text progress bar is displayed while processing MS2 groups.
#'
#' @return
#' A data frame containing one selected candidate per MS2 spectrum
#' when scoring succeeds, with groups ordered by their first
#' appearance in `matching_results`.
#'
#' Input columns are preserved. For groups containing multiple
#' candidates, a numeric `score` column records the selected
#' candidate's cosine similarity. An existing `score` column is
#' overwritten for these scored groups.
#'
#' Single-candidate groups are returned unchanged. If the input has
#' no `score` column and scored groups are also present, their output
#' scores are filled with `NA`. If every group has one candidate,
#' no new `score` column is created.
#'
#' Scores for discarded candidates, assembled chemical formulas,
#' theoretical isotope patterns, and comparison vectors are not
#' returned. No dedicated result class is assigned.
#'
#' @seealso [MatchMS2SpectraToGlycanCompositions()],
#'   [ConstructGlycanLibrary()], [enviPat::isopattern()].
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # constructed is the complete output of ConstructGlycanLibrary().
#' # diagnostic_result contains the diagnostic MS2 screening results.
#' # ms_data_raw contains the corresponding original MS1 spectra.
#'
#' # Keep multiple precursor-mass candidates for isotope comparison
#' candidates <- MatchMS2SpectraToGlycanCompositions(
#'   diagnostic_ms2_result = diagnostic_result,
#'   glycan_library = constructed,
#'   mass_error_tolerance_type = "ppm_tolerance",
#'   mass_error_tolerance_val = 10,
#'   max_possible_candidates_num = 3
#' )
#'
#' # Use built-in component formulas.
#' # Window distances below are illustrative and use m/z units.
#' selected <- ValidateGlycanCompositionByIsotopePattern(
#'   matching_results = candidates,
#'   molecular_formula = NULL,
#'   constructed_glycan_library = constructed,
#'   bin_width = 0.3,
#'   threshold_iso_probalility = 0.01,
#'   ms_data_raw = ms_data_raw,
#'   ms1_window_left = 1,
#'   ms1_window_right = 4
#' )
#'
#' head(selected)
#'
#' # A score column is created only when multi-candidate groups are scored
#' if ("score" %in% names(selected)) {
#'   selected[, c("ms2_spectrum_id", "glycan_composition", "score")]
#' }
#' }
ValidateGlycanCompositionByIsotopePattern = function(
  matching_results, 
  molecular_formula, 
  constructed_glycan_library, 
  bin_width = 0.3, 
  threshold_iso_probalility = 0.01, 
  ms_data_raw,
  ms1_window_left, ms1_window_right

) {
  reference_monosaccharides_molecular_formula = get_reference_monosaccharides_molecular_formula()
  reference_adducts_molecular_formula = get_reference_adducts_molecular_formula()
  reference_derivatizations_molecular_formula = get_reference_derivatizations_molecular_formula()

  molecular_formula_all = c(
    molecular_formula,
    reference_monosaccharides_molecular_formula,
    reference_adducts_molecular_formula,
    reference_derivatizations_molecular_formula
  )

  used_molecular_names = c(
    names(constructed_glycan_library$monosaccharides$monosaccharide_range),
    names(constructed_glycan_library$adducts$adducts_ranges),
    names(constructed_glycan_library$modifications$modification_ranges), 
    constructed_glycan_library$derivatization_type
  )

  required_cols = c("ms1_spectrum_id", "ms2_spectrum_id")

  colnames_info = colnames(matching_results)

  if (!all(required_cols %in% colnames_info)) {
    stop("Missing required columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }

  ms2_id_all = unique(matching_results[['ms2_spectrum_id']])

  ms_id_map = stats::setNames(seq_along(ms_data_raw), ms_data_raw$spectrumId)

  ms2_id_num = length(ms2_id_all)

  # progress bar
  pb = utils::txtProgressBar(min = 0, max = ms2_id_num, style = 3)
  on.exit(close(pb))

  res_list = list()
  for (i in seq_along(ms2_id_all)) {
    # progress bar
    utils::setTxtProgressBar(pb, i)

    current_ms2_id = ms2_id_all[i]

    current_matching_info = dplyr::filter(matching_results, ms2_spectrum_id == current_ms2_id)

    if (dim(current_matching_info)[1] == 1) {

      res_list[[length(res_list) + 1]] = current_matching_info
      next

    } else {

      current_ms1_id = unique(current_matching_info[['ms1_spectrum_id']])

      current_precursor_mz = unique(current_matching_info[['ms2_precursor_mz']])

      if (length(current_ms1_id) != 1) stop("ms1_spectrum_id is not unique for ms2_id=", current_ms2_id, call.=FALSE)
      if (length(current_precursor_mz) != 1) stop("ms2_precursor_mz is not unique for ms2_id=", current_ms2_id, call.=FALSE)

      ms1_window_left_range = current_precursor_mz - ms1_window_left
      ms1_window_right_range = current_precursor_mz + ms1_window_right

      this_index = which(
        as.character(ms_data_raw$spectrumId) == as.character(current_ms1_id)
      )

      current_ms1_data = as.data.frame(Spectra::peaksData(ms_data_raw)[[this_index]])

      ms1_isolation_window_data = dplyr::filter(current_ms1_data, mz >= ms1_window_left_range & mz <= ms1_window_right_range)

      molecular_list = molecular_formula_all[used_molecular_names]

      GetIsotopePatternCosineScore = function(matching_info, mol_list, threshold_iso, ms1_window_data, bin_w) {
        utils::data(isotopes, package = "enviPat", envir = environment())
        utils::data(adducts, package = "enviPat", envir = environment())

        ms1_window_tic = sum(ms1_window_data[['intensity']])
        ms1_window_data = dplyr::mutate(ms1_window_data, relative_intensity = intensity/ms1_window_tic)

        GetCenteredVector = function(mz_value, intensity_value, centers_value) {
          # attention!!!!!: mz_value, intensity_value must be in same order
          if (length(mz_value) != length(intensity_value)) {
            stop("mz_value and intensity_value must have the same length.")
          }

          center_idx_map = sapply(mz_value, function(m) {
            which.min(abs(centers_value - m))
          })

          center_vector = data.frame(mz = centers_value, intensity = rep(0, length(centers_value)))
          unique_map_idx = unique(center_idx_map)

          for (idx in unique_map_idx) {
            total_int = sum(intensity_value[center_idx_map == idx])
            center_vector$intensity[idx] = total_int
          }
          return(center_vector)
        }

        cosine_sim = function(a, b) {
          sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
        }
        cosine_value = numeric(nrow(matching_info))

        for (i in seq_len(nrow(matching_info))) {

          current_row = matching_info[i, ]
          current_charge = current_row[['charge_state']]

          molecular_formula = NULL

          for (mol_name in names(mol_list)) {
            current_formula_num = current_row[[mol_name]]

            if (current_formula_num == 0) {
              next
            }

            current_base_formula = mol_list[[mol_name]]
            current_comp_formula = enviPat::multiform(current_base_formula, current_formula_num)

            if (is.null(molecular_formula)) {
              molecular_formula = current_comp_formula
            } else {
              molecular_formula = enviPat::mergeform(molecular_formula, current_comp_formula)
            }
          }

          # isotopic distribution calculation
          molecular_formula = enviPat::check_chemform(isotopes, molecular_formula)
          new_molecular_formula = molecular_formula$new_formula

          molecular_iso_pattern = enviPat::isopattern(isotopes, new_molecular_formula,
                                                      charge = current_charge, rel_to = 2, threshold = threshold_iso, verbose = F)

          theoretical_iso_pattern = as.data.frame(molecular_iso_pattern[[1]])

          mz_min = min(theoretical_iso_pattern[['m/z']], ms1_window_data[['mz']])
          mz_max = max(theoretical_iso_pattern[['m/z']], ms1_window_data[['mz']])

          centers = seq(from = mz_min, to = mz_max, by = bin_w)

          ms1_vector = GetCenteredVector(mz_value = ms1_window_data[['mz']], intensity_value = ms1_window_data[['relative_intensity']], centers_value = centers)
          theoretical_iso_vector = GetCenteredVector(mz_value = theoretical_iso_pattern[['m/z']], intensity_value = theoretical_iso_pattern[['abundance']], centers_value = centers)

          similarity_value = cosine_sim(ms1_vector$intensity, theoretical_iso_vector$intensity)

          cosine_value[i] = similarity_value
        }

        #matching_info = dplyr::slice_max(matching_info, cosine_value, n = 1, with_ties = T)
        matching_info <- matching_info |>
          dplyr::mutate(score = cosine_value) |>
          dplyr::slice_max(score, n = 1, with_ties = FALSE)

        return(matching_info)
      }

      best_matching_info = GetIsotopePatternCosineScore(matching_info = current_matching_info,
                                                            mol_list = molecular_list,
                                                            threshold_iso = threshold_iso_probalility,
                                                            ms1_window_data = ms1_isolation_window_data,
                                                            bin_w = bin_width)

      res_list[[length(res_list) + 1]] = best_matching_info
    }
  }
  validate_spectrum_matching_info = dplyr::bind_rows(res_list)
  return(validate_spectrum_matching_info)
}











