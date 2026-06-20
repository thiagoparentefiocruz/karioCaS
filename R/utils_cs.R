# ==============================================================================
# CONFIDENCE SCORE (CS) - CANONICAL UNIT HELPERS
# ==============================================================================
#
# karioCaS stores every Confidence Score internally as an INTEGER PERCENT in the
# closed range [0, 100] (e.g. 0, 40, 90, 100). Two helpers translate the two
# places where CS enters the package into that canonical unit:
#
#   * .cs_filename_to_percent() : parses the numeric token of a `SAMPLE_CSxx`
#       filename. Accepts the legacy zero-padded "tenths" code (CS00..CS10,
#       where CS09 = 0.9 and CS10 = 1.0), an explicit percentage (CS40, CS100),
#       or an explicit decimal fraction (CS0.9, CS0.05, CS1.0).
#
#   * .cs_arg_to_percent() : parses a CS value supplied by the user as a function
#       argument (e.g. confidence_score = 1.0, CS_B = 40). Values in (0, 1] are
#       read as Kraken fractions (1.0 -> 100), values > 1 as percentages
#       (40 -> 40). This keeps existing percent-style arguments working while
#       making the natural `1.0` request for maximum stringency behave.
#
# Both return an integer in [0, 100], or NA_real_ when the input is missing,
# unparseable, or out of range.
# ==============================================================================

#' Convert a filename CS token to canonical integer percent
#' @param token Character or numeric scalar (the part after `_CS`).
#' @return Integer percent in [0, 100], or NA_real_ if invalid.
#' @noRd
.cs_filename_to_percent <- function(token) {
    token_chr <- as.character(token)
    val <- suppressWarnings(as.numeric(token_chr))
    if (length(val) != 1 || is.na(val)) {
        return(NA_real_)
    }
    has_decimal <- grepl("\\.", token_chr)
    pct <- if (has_decimal || (val > 0 && val < 1)) {
        val * 100 # decimal fraction: 0.9 -> 90, 0.05 -> 5, 1.0 -> 100
    } else if (val <= 10) {
        val * 10 # legacy tenths code: 0 -> 0, 9 -> 90, 10 -> 100
    } else {
        val # explicit percentage: 40 -> 40, 100 -> 100
    }
    pct <- round(pct)
    if (pct < 0 || pct > 100) {
        return(NA_real_)
    }
    as.integer(pct)
}

#' Convert a user-supplied CS argument to canonical integer percent
#' @param value Numeric or character scalar (a CS function argument).
#' @return Integer percent in [0, 100], or NA_real_ if invalid.
#' @noRd
.cs_arg_to_percent <- function(value) {
    val <- suppressWarnings(as.numeric(value))
    if (length(val) != 1 || is.na(val)) {
        return(NA_real_)
    }
    pct <- if (val > 0 && val <= 1) {
        val * 100 # Kraken fraction: 1.0 -> 100, 0.9 -> 90
    } else {
        val # percentage: 40 -> 40, 90 -> 90, 0 -> 0
    }
    pct <- round(pct)
    if (pct < 0 || pct > 100) {
        return(NA_real_)
    }
    as.integer(pct)
}

#' Safe maximum that tolerates empty / all-NA input
#' @param x Numeric vector.
#' @param default Value returned when x has no usable entries.
#' @return A length-one numeric.
#' @noRd
.safe_max <- function(x, default = 0) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
        return(default)
    }
    max(x)
}
