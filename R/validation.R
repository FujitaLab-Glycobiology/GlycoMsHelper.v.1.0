validate_single_positive_number = function(x, name, is_integers) {
    # length
    if (length(x) != 1) {
        stop(paste("The", name, "must be a single numeric value."))
    }
    # numeric
    if (!is.numeric(x)) {
        stop(paste("The", name, "must be numeric."))
    }
    # not NA
    if (is.na(x)) {
        stop(paste("The", name, "must not be NA."), call. = FALSE)
    }
    # finiteness
    if (!is.finite(x)) {
        stop(paste("The", name, "must be finite."))
    }
    # positive integers
    if (is_integers) {
        if (x <= 0 || x %% 1 != 0) {
            stop(paste("The", name, "must be a positive integer."))
        }
    } else {
        if (x <= 0) {
            stop(paste("The", name, "must be a positive number."))
        }
    }
}

validate_single_nonnegative_number = function(x, name, is_integers) {
    # length
    if (length(x) != 1) {
        stop(paste("The", name, "must be a single numeric value."))
    }
    # numeric
    if (!is.numeric(x)) {
        stop(paste("The", name, "must be numeric."))
    }
    # not NA
    if (is.na(x)) {
        stop(paste("The", name, "must not be NA."), call. = FALSE)
    }
    # finiteness
    if (!is.finite(x)) {
        stop(paste("The", name, "must be finite."))
    }
    # non-negative integers
    if (is_integers) {
        if (x < 0 || x %% 1 != 0) {
            stop(paste("The", name, "must be a non-negative integer."))
        }
    } else {
        if (x < 0) {
            stop(paste("The", name, "must be a non-negative number."))
        }
    }
}

validate_single_negative_number = function(x, name, is_integers) {
    # length
    if (length(x) != 1) {
        stop(paste("The", name, "must be a single numeric value."))
    }
    # numeric
    if (!is.numeric(x)) {
        stop(paste("The", name, "must be numeric."))
    }
    # not NA
    if (is.na(x)) {
        stop(paste("The", name, "must not be NA."), call. = FALSE)
    }
    # finiteness
    if (!is.finite(x)) {
        stop(paste("The", name, "must be finite."))
    }
    # negative integers
    if (is_integers) {
        if (x >= 0 || x %% 1 != 0) {
            stop(paste("The", name, "must be a negative integer."))
        }
    } else {
        if (x >= 0) {
            stop(paste("The", name, "must be a negative number."))
        }
    }
}




validate_single_number = function(x, name, is_integers) {
    # length
    if (length(x) != 1) {
        stop(paste("The", name, "must be a single numeric value."))
    }
    # numeric
    if (!is.numeric(x)) {
        stop(paste("The", name, "must be numeric."))
    }
    # not NA
    if (is.na(x)) {
        stop(paste("The", name, "must not be NA."), call. = FALSE)
    }
    # finiteness
    if (!is.finite(x)) {
        stop(paste("The", name, "must be finite."))
    }
    # integers
    if (is_integers) {
        if (x %% 1 != 0) {
            stop(paste("The", name, "must be an integer."))
        }
    }
}








validate_positive_number_vector <- function(x, name, is_integers, is_named_vector) {
    # not empty
    if (length(x) == 0) {
        stop(paste("The", name, "must not be empty."))
    }
    # numeric
    if (!is.numeric(x)) {
        stop(paste("The", name, "must be numeric."))
    }
    # finiteness
    if (!all(is.finite(x))) {
        stop(paste("All values in", name, "must be finite."))
    }
    # positive 
    if (is_integers) {
        # integers
        if (any(x <= 0) || any(x %% 1 != 0)) {
            stop(paste("All values in", name, "must be positive integers."))
        }
    } else {
        # non-integers
        if (any(x <= 0)) {
            stop(paste("All values in", name, "must be positive."))
        }
    }
    # no NA
    if (any(is.na(x))) {
        stop(paste("The", name, "must not contain NA values."))
    }
    # named vector
    if (is_named_vector) {
        # names are required
        if (is.null(names(x)) || any(names(x) == "")) {
            stop(paste("The", name, "must be a named vector."))
        }
        # all names must be unique
        if (any(duplicated(names(x)))) {
            stop(paste("All names in", name, "must be unique."))
        }
    }
}

validate_nonnegative_number_vector <- function(x, name, is_integers, is_named_vector) {
    # not empty
    if (length(x) == 0) {
        stop(paste("The", name, "must not be empty."))
    }
    # numeric
    if (!is.numeric(x)) {
        stop(paste("The", name, "must be numeric."))
    }
    # finiteness
    if (!all(is.finite(x))) {
        stop(paste("All values in", name, "must be finite."))
    }
    # non-negative
    if (is_integers) {
        #  integers
        if (any(x < 0) || any(x %% 1 != 0)) {
            stop(paste("All values in", name, "must be non-negative integers."))
        }
    } else {
        # non-integers
        if (any(x < 0)) {
            stop(paste("All values in", name, "must be non-negative."))
        }
    }
    # no NA
    if (any(is.na(x))) {
        stop(paste("The", name, "must not contain NA values."))
    }

    if (is_named_vector) {
        # names are required
        if (is.null(names(x)) || any(names(x) == "")) {
            stop(paste("The", name, "must be a named vector."))
        }
        # all names must be unique
        if (any(duplicated(names(x)))) {
            stop(paste("All names in", name, "must be unique."))
        }
    }
}





validate_single_logical <- function(x, argument_name) {
    if (!is.logical(x) || length(x) != 1 || is.na(x)) {
        stop(paste0("`", argument_name, "` must be a single TRUE or FALSE value.")) 
    }
}




