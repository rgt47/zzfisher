# syntax=docker/dockerfile:1.4
# zzcollab Dockerfile v2.4.0

ARG BASE_IMAGE=rocker/tidyverse
ARG R_VERSION=4.4.2
ARG USERNAME=analyst

FROM ${BASE_IMAGE}:${R_VERSION}

ARG USERNAME=analyst
ARG DEBIAN_FRONTEND=noninteractive

# RENV_CONFIG_REPOS_OVERRIDE forces renv to use Posit PPM binaries
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 TZ=UTC \
    RENV_PATHS_CACHE=/home/${USERNAME}/.cache/R/renv \
    RENV_CONFIG_REPOS_OVERRIDE="https://packagemanager.posit.co/cran/__linux__/noble/latest" \
    ZZCOLLAB_CONTAINER=true

# No additional system dependencies required

# Configure R to use Posit Package Manager for pre-compiled binaries
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))' \
        >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))' \
        >> /usr/local/lib/R/etc/Rprofile.site

# Install renv and restore packages from lockfile (using PPM binaries)
RUN R -e "install.packages('renv')"
RUN mkdir -p /home/${USERNAME}/.cache/R/renv && chmod 777 /home/${USERNAME}/.cache/R/renv
COPY renv.lock renv.lock
RUN R -e "renv::restore()"

# Install languageserver for IDE support and yaml for R Markdown dependencies
RUN R -e "install.packages(c('languageserver', 'yaml'))"

# Create non-root user
RUN useradd --create-home --shell /bin/bash ${USERNAME} && \
    chown -R ${USERNAME}:${USERNAME} /usr/local/lib/R/site-library

USER ${USERNAME}
WORKDIR /home/${USERNAME}/project

CMD ["R", "--quiet"]
