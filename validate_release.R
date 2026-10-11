repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(repo, "README.md")))

r_files <- list.files(repo, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
stan_files <- list.files(file.path(repo, "models"), pattern = "\\.stan$",
                        full.names = TRUE)
if (!length(r_files) || !length(stan_files)) stop("R or Stan sources are missing.")

parse_errors <- lapply(r_files, function(path) {
  tryCatch({ parse(path); NULL }, error = function(e) paste(path, conditionMessage(e)))
})
parse_errors <- unlist(parse_errors)
if (length(parse_errors)) stop(paste(parse_errors, collapse = "\n"))

text_files <- c(
  list.files(repo, pattern = "\\.(R|stan|md|cff|txt)$", recursive = TRUE,
             full.names = TRUE),
  file.path(repo, ".gitignore"), file.path(repo, ".gitattributes")
)
text_files <- setdiff(normalizePath(text_files, winslash = "/", mustWork = TRUE),
                      normalizePath(file.path(repo, "validate_release.R"),
                                    winslash = "/", mustWork = TRUE))
contents <- paste(vapply(text_files, function(x) paste(readLines(x, warn = FALSE),
                                                       collapse = "\n"), character(1)),
                  collapse = "\n")
forbidden <- c("LAB-Mobile", "Rambam", "NITSAN", "SubID", "[A-Z]:[/\\\\]Users")
hits <- forbidden[vapply(forbidden, grepl, logical(1), x = contents,
                         ignore.case = FALSE, perl = TRUE)]
if (length(hits)) stop("Private/internal token(s) detected: ", paste(hits, collapse = ", "))

private_inputs <- list.files(file.path(repo, "data"), pattern = "\\.(csv|tsv|xlsx)$",
                             full.names = TRUE)
if (length(private_inputs)) stop("Participant-level input files must not be committed.")

cat("Release validation passed:\n")
cat(" -", length(r_files), "R files parsed\n")
cat(" -", length(stan_files), "Stan model files present\n")
cat(" - no forbidden local/path/source-ID tokens found\n")
cat(" - no participant-level data files present\n")
