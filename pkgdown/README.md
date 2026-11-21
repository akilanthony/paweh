# Website maintenance

The website uses Bootstrap 5 through pkgdown. Edit the homepage in `index.md`
and the stylesheet in `extra.css`. Navigation and topic groups are defined in
the repository’s `_pkgdown.yml`.

For a clean preview, work from a temporary copy of the package source that
includes `pkgdown/`, `man/figures/`, `_pkgdown.yml` and all eight vignettes.
Exclude `AGENTS.md` and `CRAN_PREPARATION_AUDIT.md` from the copy because
pkgdown renders top-level Markdown pages. With pkgdown and the package
dependencies installed, run this from the temporary source directory.

```r
pkgdown::build_site_github_pages(install = TRUE, new_process = FALSE)
```

Open `docs/index.html`. To test site search or other features that request
local files, serve `docs/` over HTTP on localhost. For example, run
`python3 -m http.server 8000 --bind 127.0.0.1 --directory docs` and open
`http://127.0.0.1:8000/`.

## Publishing after review

Pushes to `develop` and pull requests build the site. Only a manual
`workflow_dispatch` on `refs/heads/develop` can run the publishing job.
It uses the site artifact from that same run and appends a commit to
`gh-pages` with force-push disabled.

GitHub requires the [dispatch workflow on the default branch](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow).
The default branch is currently `main`, where the file is absent. Register
the gated workflow there as part of a separate approved publishing step,
then select `develop` in the manual dispatch UI. Keep the Pages source set
to the root of the `gh-pages` branch.

Configure the repository Actions secret `PAWEH_PAGES_TOKEN` with a credential
limited to this repository and permitted to write to `gh-pages`. Branch-based
GitHub Pages builds require a push with a credential other than the workflow’s
default `GITHUB_TOKEN`. See the
[GitHub Pages publishing-source documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).
