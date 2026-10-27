get_reference_monosaccharides <- function() {
    c(
        Hex = 162.0528,
        HexNAc = 203.0794,
        dHex = 146.0579,
        Neu5Ac = 291.0954,
        HexA = 176.03209,
        Neu5Gc = 307.0903,
        Pentose = 132.0423,
        KDN = 250.0689, 
        EtNP = 123.0085291, 
        AHM = 144.0422587
    )
}

get_reference_adducts <- function() {
    data.frame(
        adduct = c("H", "Na", "K", "Li", "Mg"),
        mass = c(1.00727, 22.98922070, 38.963707, 7.016003, 23.985042),
        charge = c(1, 1, 1, 1, 2),
        stringsAsFactors = FALSE
    )
}

get_reference_derivatizations <- function() {
    c(
        ProA = 237.1841124,
        AB = 138.0793129,
        PA = 96.068761
    )
}

get_supported_glycan_type = function() {
    c("N_glycan", "O_glycan", "GPI", "GSL", "HMO")
}


get_glycan_type_std_monosaccharides <- function(glycan_type) {
    if (glycan_type == "N_glycan") {
        return(c('Hex', 'HexNAc', 'dHex', 'Neu5Ac', 'HexA', 'Neu5Gc'))
    } else if (glycan_type == "O_glycan") {
        return(c('Hex', 'HexNAc', 'dHex', 'Neu5Ac', 'HexA', 'Neu5Gc'))
    } else if (glycan_type == "GPI") {
        return(c('Hex', 'HexNAc', 'Neu5Ac', 'Neu5Gc', 'EtNP', 'AHM'))
    } else if (glycan_type == "GSL") {
        return(c('Hex', 'HexNAc', 'dHex', 'Neu5Ac', 'Neu5Gc'))
    } else if (glycan_type == "HMO") {
        return(c('Hex', 'HexNAc', 'dHex', 'Neu5Ac'))
    } else {
        stop("Unsupported glycan type")
    }
}


get_most_abundant_iso_name <- function() {
  c(
    "1H","4He","7Li","9Be","11B","12C","14N","16O","19F","20Ne", "23Na", 
    "24Mg","27Al","28Si","31P","32S","35Cl","40Ar","39K","40Ca",
    "45Sc","48Ti","51V","52Cr","55Mn","56Fe","59Co","58Ni","63Cu","64Zn",
    "69Ga","74Ge","75As","80Se","79Br","84Kr","85Rb","88Sr","89Y","90Zr",
    "93Nb","98Mo","102Ru","103Rh","106Pd","107Ag","114Cd",
    "115In","120Sn","121Sb","130Te","127I","132Xe","133Cs","138Ba",
    "139La","140Ce","141Pr","142Nd","152Sm","153Eu","158Gd",
    "159Tb","164Dy","165Ho","166Er","169Tm","174Yb","175Lu","180Hf",
    "181Ta","184W","187Re","192Os","193Ir","195Pt","197Au","202Hg",
    "205Tl","208Pb","209Bi", "232Th", "231Pa", "238U"
  )
}


get_other_iso_name <- function() {
  c(
    "2H", "3He", "6Li", "10B", "13C", "15N", "17O", "18O", "21Ne", "22Ne", "25Mg", "26Mg", "29Si", "30Si", "33S", "34S", "35S", "36S",
    "36Cl",  "37Cl",  "36Ar",  "38Ar",  "40K" ,  "41K", "41Ca",  "42Ca", "43Ca",  "44Ca", "45Ca",  "46Ca",  "47Ca",  "48Ca",  "46Ti",  "47Ti",  "49Ti",  "50Ti",
    "50V",   "50Cr" , "53Cr",  "54Cr",  "54Fe" , "55Fe" , "57Fe",  "58Fe",  "60Ni",  "61Ni",  "62Ni",  "64Ni",  "65Cu",  "66Zn",  "67Zn",  "68Zn",  "70Zn",  "71Ga" ,
    "70Ge",  "72Ge"  ,"73Ge",  "76Ge",  "74Se" , "76Se" , "77Se",  "78Se",  "82Se",  "80Br",  "81Br",  "78Kr",  "80Kr",  "82Kr",  "83Kr",  "86Kr",  "87Rb",  "84Sr" ,
    "86Sr" , "87Sr" , "91Zr",  "92Zr",  "94Zr",  "96Zr" , "92Mo",  "94Mo",  "95Mo",  "96Mo" , "97Mo",  "100Mo", "96Ru",  "98Ru",  "99Ru",  "100Ru", "101Ru", "104Ru",
    "102Pd", "104Pd", "105Pd", "108Pd", "110Pd" ,"109Ag", "106Cd", "108Cd", "110Cd", "111Cd" ,"112Cd", "113Cd", "116Cd", "113In", "112Sn", "114Sn", "115Sn", "116Sn",
    "117Sn", "118Sn", "119Sn", "122Sn", "124Sn" ,"123Sb", "120Te", "122Te", "123Te", "124Te" ,"125Te", "126Te", "128Te", "124Xe", "126Xe", "128Xe", "129Xe", "130Xe",
    "131Xe", "134Xe", "136Xe", "130Ba", "132Ba" ,"134Ba", "135Ba", "136Ba", "137Ba", "138La" ,"136Ce", "138Ce", "142Ce", "143Nd", "144Nd", "145Nd", "146Nd", "148Nd",
    "150Nd", "144Sm", "147Sm", "148Sm", "149Sm" ,"150Sm", "154Sm", "151Eu", "152Gd", "154Gd", "155Gd", "156Gd", "157Gd", "160Gd", "156Dy", "158Dy", "160Dy", "161Dy",
    "162Dy", "163Dy", "162Er", "164Er", "167Er" ,"168Er", "170Er", "168Yb", "170Yb", "171Yb", "172Yb", "173Yb", "176Yb", "176Lu", "174Hf", "176Hf", "177Hf", "178Hf",
    "179Hf", "180Ta", "180W",  "182W",  "183W"  ,"186W",  "185Re", "184Os", "186Os", "187Os", "188Os", "189Os", "190Os", "191Ir", "190Pt", "192Pt", "194Pt", "196Pt",
    "198Pt", "196Hg", "198Hg", "199Hg", "200Hg" ,"201Hg", "204Hg", "203Tl", "204Pb", "206Pb", "207Pb", "234U" , "235U"
  )
}



get_denoising_methods <- function() {
    c("adaptive", "spline", "segmented", "quantile", "fixed")
}

get_transform_methods <- function() {
    c("log2", "sqrt", "asinh", "none")
    # c("log2", "log10", "log", "sqrt", "none")
}



get_mass_error_tolerance_types <- function() {
    c("ppm_tolerance", "absolute_tolerance")
}


get_ms_data_filtering_methods <- function() {
    c("mean_sd", "quantile", "start_end")
}


get_supported_ms_data_filtering_variables <- function() {
    c("peaksCount", "totIonCurrent", "rtime")
}






get_reference_monosaccharides_molecular_formula <- function() {
    c(
        Hex = 'C6H10O5',
        HexNAc = 'C8H13N1O5',
        dHex = 'C6H10O4',
        Neu5Ac = 'C11H17N1O8',
        HexA = 'C6H8O6',
        Neu5Gc = 'C11H17N1O9',
        Pentose = 'C5H8O4',
        KDN = 'C9H14O8',
        EtNP = 'C2H6N1O3P1',
        AHM = 'C6H8O4'
    )
}

get_reference_adducts_molecular_formula <- function() {
    c(
        H = 'H1', 
        Na = 'Na1',
        K = 'K1',
        Li = 'Li1',
        Mg = 'Mg1'
    )
}

get_reference_derivatizations_molecular_formula <- function() {
    c(
        ProA = 'C13H23N3O1',
        AB = 'C7H10N2',
        PA = 'C5H8N2'
    )
}