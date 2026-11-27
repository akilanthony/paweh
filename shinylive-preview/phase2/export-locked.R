# Run from the repository root. CACHE can be prepopulated for offline repeat builds.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop("Supply ARCHIVE DESTINATION CACHE [PROVENANCE_JSON]", call. = FALSE)
}
archive <- normalizePath(args[[1L]], mustWork = TRUE)
destination <- args[[2L]]
cache <- args[[3L]]
lock <- jsonlite::read_json("shinylive-preview/phase2/lock.json")
stopifnot(packageVersion("shinylive") == lock$shinylive)
stopifnot("paweh/DESCRIPTION" %in% utils::untar(archive, list = TRUE, tar = "internal"))

sum_file <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  unclass(as.character(openssl::sha256(con)))
}

expected_archive_hash <- lock$portable_paweh_sha256
if (length(args) == 4L) {
  provenance <- jsonlite::read_json(args[[4L]])
  stopifnot(identical(provenance$paweh_source, lock$source_commit),
            identical(provenance$wasm_r, lock$wasm_r))
  expected_archive_hash <- provenance$wrapper_archive_sha256
}
stopifnot(identical(sum_file(archive), expected_archive_hash))

verified_input <- function(url, sha, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(path)) {
    status <- system2("curl", c("--fail", "--location", "--silent", "--show-error",
      "--retry", "2", "--output", shQuote(path), shQuote(url)))
    if (status != 0) stop("Download failed: ", url, call. = FALSE)
  }
  if (!identical(sum_file(path), sha)) stop("Checksum mismatch: ", path, call. = FALSE)
  path
}

Sys.setenv(R_USER_CACHE_DIR = tempfile("shinylive-verified-assets-"))
asset_dir <- rappdirs::user_cache_dir("shinylive")
asset_bundle <- verified_input(lock$web_assets_url, lock$web_assets_sha256,
  file.path(cache, basename(lock$web_assets_url)))
dir.create(asset_dir, recursive = TRUE)
stopifnot(utils::untar(asset_bundle, exdir = asset_dir) == 0)

appdir <- tempfile("paweh-shinylive-")
dir.create(appdir)
stopifnot(file.copy("shinylive-preview/app.R", appdir))
stopifnot(file.copy(archive, file.path(appdir, "paweh-library.tar.gz")))
shinylive::export(appdir, destination, wasm_packages = FALSE,
  assets_version = lock$web_assets, template_params = list(title = "paweh"))

webr <- file.path(destination, "shinylive", "webr")
records <- lapply(lock$packages, function(p) {
  file <- verified_input(p$url, p$sha256, file.path(cache, basename(p$path)))
  installed <- file.path(webr, p$path)
  dir.create(dirname(installed), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file, installed))
  list(name = p$name, version = p$version, cached = TRUE,
       assets = list(list(filename = basename(p$path), url = p$url, version = p$version)),
       type = "package", ref = p$ref, path = p$path)
})
names(records) <- vapply(records, `[[`, "", "name")
saveRDS(records, file.path(webr, "packages", "metadata.rds"))
