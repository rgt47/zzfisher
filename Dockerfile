# syntax=docker/dockerfile:1.4
# zzcollab Dockerfile v0.1.0

# BASE_IMAGE is parsed out of this file by the project Makefile ('make r'
# derives the profile label from it); keep it even though the FROM below uses
# a fully-substituted literal and does not reference the ARG.
ARG BASE_IMAGE=rocker/verse

FROM rocker/verse:4.6.0@sha256:0b50fdee9288723b5a6802502341f0429ab05edb9db04d57958f49e18d3ea883

# OCI image labels for reproducibility provenance and tooling integration.
# base_digest records the resolved sha256 of the rocker base at build time;
# ppm_snapshot records the dated PPM URL used to pin package binaries.
LABEL org.opencontainers.image.created="2026-07-06T00:28:07Z" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      zzcollab.template.version="0.1.0" \
      zzcollab.r.version="4.6.0" \
      zzcollab.base.image="rocker/verse:4.6.0" \
      zzcollab.base.digest="sha256:0b50fdee9288723b5a6802502341f0429ab05edb9db04d57958f49e18d3ea883" \
      zzcollab.ppm.snapshot="2026-07-05" \
      zzcollab.install.mode="renv"

ARG USERNAME=analyst
ARG DEBIAN_FRONTEND=noninteractive

# RENV_PATHS_LIBRARY is outside the project bind-mount so the baked library
# is not shadowed at runtime. ZZCOLLAB_AUTO_RESTORE=false disables the
# startup restore so the image library is authoritative.
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 TZ=UTC \
    RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_PATHS_CACHE=/opt/renv/cache \
    RENV_CONFIG_REPOS_OVERRIDE="https://packagemanager.posit.co/cran/__linux__/noble/2026-07-05" \
    ZZCOLLAB_CONTAINER=true \
    ZZCOLLAB_INSTALL_MODE=renv \
    ZZCOLLAB_AUTO_RESTORE=false

# No additional system dependencies required

# Configure R to use Posit Package Manager for pre-compiled binaries
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/2026-07-05"))' \
        >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))' \
        >> /usr/local/lib/R/etc/Rprofile.site

# Install R dev tooling only (languageserver/styler/lintr, config-gated;
# never in renv.lock). Project packages come from renv::restore(), not here.
RUN R -e "install.packages(c('languageserver'), Ncpus = max(1L, parallel::detectCores()))"

# Pre-bake the LaTeX package closure at build time (as root, where tlmgr can
# write) so PDF rendering works for the non-root user with no runtime install.
# The closure is installed in one bulk tlmgr pass: relying on tinytex to
# discover packages lazily during a render installs them one at a time, and
# each missing .sty triggers a full pdflatex recompile (~14s x ~24 packages,
# ~5 min of build). A single tlmgr call fetches them in ~20-30s; the render
# below then finds everything present and acts as a fast smoke test that also
# self-heals any package this list omits.
RUN <<'WARMUP'
set -eu
R -e "tinytex::tlmgr_install(c('amsfonts','booktabs','setspace','multirow','wrapfig','float','colortbl','pdflscape','tabu','varwidth','threeparttable','threeparttablex','environ','trimspaces','ulem','makecell','mathtools','fancyhdr','caption','enumitem','fp','pgf','pgfplots','siunitx','lineno'))"
d=/tmp/texwarmup
mkdir -p "$d"
cat > "$d/01-report.Rmd" <<'RMD'
---
title: warm-up report
output:
  bookdown::pdf_document2:
    number_sections: true
header-includes:
  - \usepackage{setspace}
---
# Section
Math $\alpha \in \mathbb{R}$ and $\sum_{i=1}^{n} x_i^2$.
```{r}
knitr::kable(head(mtcars, 3), booktabs = TRUE)
```
RMD
cat > "$d/02-kitchensink.Rmd" <<'RMD'
---
title: warm-up kitchen sink
output:
  pdf_document:
    latex_engine: xelatex
    extra_dependencies:
      - booktabs
      - longtable
      - array
      - multirow
      - wrapfig
      - float
      - colortbl
      - pdflscape
      - tabu
      - threeparttable
      - threeparttablex
      - ulem
      - makecell
      - xcolor
      - amsmath
      - amssymb
      - amsfonts
      - mathtools
      - hyperref
      - geometry
      - fancyhdr
      - caption
      - subcaption
      - graphicx
      - multicol
      - setspace
      - enumitem
      - tikz
      - pgfplots
      - siunitx
---
# Kitchen sink
Math $\mathcal{N}(\mu, \sigma^2)$.
RMD
R -e 'options(tinytex.install_packages = TRUE); for (f in list.files("/tmp/texwarmup", pattern = "[.]Rmd$", full.names = TRUE)) rmarkdown::render(f, quiet = TRUE)'
rm -rf "$d"
WARMUP

# Dependency install (self-adapting, INSTALL_MODE=renv). The block
# below is emitted by generation-time branch on renv.lock presence. In renv
# mode, tools_install above runs BEFORE renv::init so IDE tools are in the
# system library; renv::init then activates renv and routes later installs to
# RENV_PATHS_LIBRARY.
RUN R -e "install.packages('renv')"
# 0777 so the non-root run user can hydrate/snapshot into the library (F-2);
# single-user research container, so world-writable here is acceptable.
RUN mkdir -p /opt/renv/library /opt/renv/cache && chmod 777 /opt/renv/library /opt/renv/cache
COPY renv.lock renv.lock
# RENV_LOCK_HASH is passed by the builder as a digest of renv.lock. Declaring
# it here and referencing it in the RUN below makes the restore layer's cache
# key depend on the lockfile content, so renv::restore() re-runs whenever
# renv.lock changes. This guards against BuildKit serving a stale restore
# layer, which would otherwise bake a library that silently diverges from the
# lockfile (and from the image's content-addressable hash label).
ARG RENV_LOCK_HASH=unknown
# renv::init creates the platform-specific library directory structure that
# renv::restore() requires to link packages from the cache.
RUN echo "renv.lock hash: ${RENV_LOCK_HASH}" &&     R -e "renv::init(bare=TRUE, force=TRUE, restart=FALSE); renv::restore(exclude = 'renv')"

# Install zzrenvcheck as a validation tool (system library, outside project renv).
# Installed post-build via make install-zzrenvcheck to avoid GitHub/network
# issues during docker build on cloud-mounted filesystems.


# Create non-root user, in the 'staff' group. rocker/verse owns its TeX tree
# (/opt/texlive, /usr/local/texlive) as root:staff and makes it group-writable,
# so a render that installs LaTeX packages at run time (tinytex) needs the run
# user to be in 'staff'; otherwise tlmgr/fmtutil fail with permission errors.
# Own the renv library AND cache (populated as root by the restore above) so the
# run user can hydrate/snapshot into them; the earlier chmod is non-recursive
# and predates the restore, so it does not cover the package subdirectories (F-2).
RUN useradd --create-home --shell /bin/bash --groups staff ${USERNAME} && \
    chown -R ${USERNAME}:${USERNAME} /usr/local/lib/R/site-library /opt/renv

USER ${USERNAME}
WORKDIR /home/${USERNAME}/project

CMD ["R", "--quiet"]
