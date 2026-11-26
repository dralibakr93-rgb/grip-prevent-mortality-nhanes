# Checks the raw files in data/raw/ against the SHA-256 checksums in data/raw_data_manifest.txt.
# With --download, a missing or mismatching file is downloaded again from NCHS.
# Usage, from the repository root:
#   Rscript data/download_raw_data.R              # check only
#   Rscript data/download_raw_data.R --download   # download missing files
# NCHS occasionally reissues files; a file whose checksum differs from the manifest is reported,
# because the results were produced with the listed versions.

download <- "--download" %in% commandArgs(trailingOnly = TRUE)
man <- readLines("data/raw_data_manifest.txt")
man <- man[!grepl("^#", man) & nzchar(man)]
parts <- strsplit(man, " | ", fixed = TRUE)
ok <- TRUE
for (p in parts) {
  path <- file.path("data", "raw", p[1])
  sha <- p[3]
  url <- sub(" .*$", "", p[4])
  have <- file.exists(path) && identical(digest::digest(file = path, algo = "sha256"), sha)
  if (!have && download && grepl("^https://", url)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    message("Downloading ", url)
    utils::download.file(url, path, mode = "wb", quiet = TRUE)
    have <- identical(digest::digest(file = path, algo = "sha256"), sha)
    if (!have) message("  checksum differs from stored version: ", path)
  }
  message(sprintf("%-60s %s", p[1], if (have) "OK" else "MISSING/MISMATCH"))
  ok <- ok && have
}
if (!ok) quit(status = 1)
message("All raw inputs verified.")
