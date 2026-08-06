# ==== Design principles: see R/core.R ====
# Internal-only helpers (never exported). Prefixed with `.` by convention
# so users never see them in autocomplete, and no input validation is
# duplicated here -- callers in core.R already validated their inputs.

.utd_check_names <- function(x, what = "argument") {
  if (is.null(names(x)) || any(names(x) == "")) {
    stop(sprintf("`%s` must be fully named.", what), call. = FALSE)
  }
  invisible(TRUE)
}
