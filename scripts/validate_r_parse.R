#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  args <- sort(list.files(
    "scripts",
    pattern = "\\.[Rr]$",
    recursive = TRUE,
    full.names = TRUE
  ))
}
if (length(args) == 0L) {
  stop("No R scripts found for parse validation")
}

for (path in args) {
  parse(file = path)
}

cat(sprintf("R_PARSE_PASS %d\n", length(args)))
