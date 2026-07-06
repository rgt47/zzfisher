# zzcollab Makefile v0.1.0
# Docker-first workflow for reproducible research

# Fail-fast recipes: a command failing (or an unset variable) aborts the
# recipe instead of silently continuing. pipefail is intentionally omitted
# because the `r` target uses `grep | head` pipelines that exit via SIGPIPE.
SHELL := /bin/bash
.SHELLFLAGS := -eu -c
.DEFAULT_GOAL := help

# Auto-detect from project (no manual configuration needed)
PACKAGE_NAME := $(shell basename $(CURDIR))
PROJECT_NAME := $(PACKAGE_NAME)
R_VERSION := $(shell grep 'R_VERSION=' Dockerfile 2>/dev/null | head -1 | sed 's/.*=//' || echo "4.4.0")

# Git-based versioning for reproducibility (use git SHA or date)
GIT_SHA := $(shell git rev-parse --short HEAD 2>/dev/null || echo "$(shell date +%Y%m%d)")
IMAGE_TAG = $(GIT_SHA)

# Report rendered by docker-render (override: make docker-render REPORT=path)
REPORT ?= analysis/report/report.Rmd

# Docker Hub account for team image push/pull. Override on the command line
# (make docker-push-team DOCKERHUB_ACCOUNT=youraccount) or export it.
DOCKERHUB_ACCOUNT ?=

# Raw data directory hashed by hash-data/verify-data (override if the
# project's data lives elsewhere).
RAW_DATA_DIR ?= analysis/data/raw_data

# Container execution. The mount and work path follow the Dockerfile's
# ARG USERNAME (default analyst), matching what the `r` target resolves at
# runtime, so a customised container user is honoured rather than hardcoded.
# PLATFORM is pinned so every batch container run matches CI and the cached
# test target rather than the host architecture (the real prior inconsistency).
PLATFORM       ?= linux/amd64
# Container runtime: docker (default), podman (drop-in for `run`), or apptainer
# (HPC; executes a SIF built from the project image, no --platform). The default
# below is the configured zzcollab runtime; override per-invocation with
# `make CONTAINER_RUNTIME=podman r`.
CONTAINER_RUNTIME ?= docker
SIF ?= env.sif
# Backend-aware prefix for the host-style R targets (document, build, check,
# test, ...). With a Nix flake present the environment is entered via
# `nix develop -c`; otherwise the commands run on the host (renv) directly.
# This is the same presence-branch idiom the Dockerfile and workflows use.
ifneq ($(wildcard flake.nix),)
  R_RUN := nix develop -c
else
  R_RUN :=
endif
CONTAINER_USER := $(shell grep '^ARG USERNAME=' Dockerfile 2>/dev/null | head -1 | cut -d= -f2)
CONTAINER_HOME := /home/$(if $(CONTAINER_USER),$(CONTAINER_USER),analyst)
# Apptainer execs a SIF (auto-binds the CWD; no --platform); docker/podman run
# the OCI image with explicit mounts. Both expose the project at the CWD.
ifeq ($(CONTAINER_RUNTIME),apptainer)
  DOCKER_RUN = apptainer exec $(SIF)
else
  DOCKER_RUN = $(CONTAINER_RUNTIME) run --rm --platform $(PLATFORM) \
    -v "$(CURDIR)":$(CONTAINER_HOME)/project \
    -w $(CONTAINER_HOME)/project $(PACKAGE_NAME)
endif

# Self-documenting help. Targets are described with a trailing `## text` and
# grouped by `##@ Section` markers, so the listing cannot drift from reality.
help: ## Show this help
	@awk 'BEGIN {FS = ":.*## "} \
	  /^##@/ {printf "\n%s\n", substr($$0, 5); next} \
	  /^[a-zA-Z0-9_.-]+:.*## / {printf "  %-22s %s\n", $$1, $$2}' \
	  $(MAKEFILE_LIST)

##@ Validation & formatting (deps checked on host; snapshot/style/lint use Docker)

# Dependency validation runs zzrenvcheck on the HOST, not in the container: it
# is static analysis (it reads R/, analysis/, scripts/ - and tests/, vignettes/
# in strict mode - plus DESCRIPTION and renv.lock, and queries the repos for
# versions) and never needs the container's installed library. Only
# renv::snapshot, which does need that library, runs in the container. Host
# repos are pinned to the image's PPM snapshot (parsed from the Dockerfile) so
# any version zzrenvcheck resolves matches the container; falls back to CRAN
# when no Dockerfile is present.
PPM_REPO := $(shell sed -n 's/.*RENV_CONFIG_REPOS_OVERRIDE="\([^"]*\)".*/\1/p' Dockerfile 2>/dev/null | head -1)
PPM_REPO := $(if $(PPM_REPO),$(PPM_REPO),https://cloud.r-project.org)
HOST_CHECK := Rscript -e "options(repos = c(CRAN = '$(PPM_REPO)'))" -e

ensure-zzrenvcheck: ## Install zzrenvcheck on the host if missing (validation tool)
	@Rscript -e "if (!requireNamespace('zzrenvcheck', quietly = TRUE)) { install.packages('remotes'); remotes::install_github('rgt47/zzrenvcheck', upgrade = 'never') }"

check-renv: ensure-zzrenvcheck ## Validate deps on host: strict + auto-fix (run before commit)
	$(HOST_CHECK) "zzrenvcheck::check_packages(auto_fix = TRUE, strict = TRUE)"

check-renv-no-fix: ensure-zzrenvcheck ## Validate deps on host: report only, no modifications
	$(HOST_CHECK) "zzrenvcheck::check_packages(auto_fix = FALSE, strict = TRUE)"

check-renv-no-strict: ensure-zzrenvcheck ## Validate deps on host: skip tests/ and vignettes/
	$(HOST_CHECK) "zzrenvcheck::check_packages(auto_fix = TRUE, strict = FALSE)"

snapshot: ensure-zzrenvcheck ## renv snapshot (container) then validate deps (host); run before commit
	$(DOCKER_RUN) Rscript -e "renv::hydrate(); renv::snapshot(prompt = FALSE)"
	$(HOST_CHECK) "zzrenvcheck::check_packages(auto_fix = TRUE, strict = TRUE)"

reset-lock: ensure-zzrenvcheck ## Rebuild renv.lock from a fresh host code scan, refreshing ALL versions
	$(HOST_CHECK) "zzrenvcheck::check_packages(fresh = TRUE, strict = TRUE)"

style: ## Format R code with styler (in container; needs code-quality enabled)
	$(DOCKER_RUN) Rscript -e "styler::style_dir('R'); if (dir.exists('analysis')) styler::style_dir('analysis')"

lint: ## Lint R code with lintr (in container; needs code-quality enabled)
	$(DOCKER_RUN) Rscript -e "print(lintr::lint_dir('R')); if (dir.exists('analysis')) print(lintr::lint_dir('analysis'))"

##@ Main workflow (recommended)

# Terminal: Interactive bash for vim editing.
# No check-renv prerequisite: that spins up a full container validation before
# you get a shell, and a second validation already runs after the session.
r: ## Interactive R session (Docker/Podman, nix, apptainer; host R fallback if no container)
	@if [ "$(CONTAINER_RUNTIME)" = apptainer ]; then \
		if [ ! -f $(SIF) ]; then \
			echo "No $(SIF) found - build it first: make sif"; exit 1; \
		fi; \
		echo "📦 Apptainer: entering 'apptainer shell $(SIF)'..."; \
		apptainer shell $(SIF); \
		exit 0; \
	fi; \
	if [ -f flake.nix ]; then \
		if command -v nix >/dev/null 2>&1; then \
			echo "❄️  Nix backend detected; entering 'nix develop'..."; \
			nix develop; \
			exit 0; \
		fi; \
		echo "⚠️  flake.nix present but 'nix' is not installed; skipping Nix."; \
	fi; \
	if [ ! -f Dockerfile ] || ! $(CONTAINER_RUNTIME) image inspect $(PACKAGE_NAME) >/dev/null 2>&1; then \
		echo ""; \
		if [ ! -f Dockerfile ]; then \
			echo "⚠️  No Dockerfile found - workspace not initialized."; \
		elif ! command -v $(CONTAINER_RUNTIME) >/dev/null 2>&1; then \
			echo "⚠️  '$(CONTAINER_RUNTIME)' is not installed."; \
		elif ! $(CONTAINER_RUNTIME) info >/dev/null 2>&1; then \
			echo "⚠️  '$(CONTAINER_RUNTIME)' is installed but not responding"; \
			echo "    (daemon stopped, or 'podman machine' not started?)."; \
		else \
			echo "⚠️  Container image '$(PACKAGE_NAME)' not built."; \
		fi; \
		if command -v R >/dev/null 2>&1; then \
			echo ""; \
			echo "Host R available: $$(command -v R)"; \
			echo "WARNING: host packages may differ from the container;"; \
			echo "         a host session is not guaranteed reproducible."; \
			printf "Run host R anyway? [y/N] "; \
			read _ans; \
			case "$$_ans" in \
				[Yy]*) echo ""; R; exit 0 ;; \
				*) echo "Aborted."; exit 1 ;; \
			esac; \
		fi; \
		echo ""; \
		echo "❌ No container available and no host R found."; \
		echo ""; \
		echo "Run zzcollab to create a Docker environment:"; \
		echo ""; \
		echo "  zzcollab docker                            # default profile"; \
		echo "  zzcollab docker --profile analysis         # tidyverse"; \
		echo "  zzcollab docker --base-image rocker/verse  # LaTeX/Quarto via the rocker/verse base image"; \
		echo ""; \
		echo "See: zzcollab docker --help for all options"; \
		echo ""; \
		exit 1; \
	fi; \
	echo "🔍 Checking workspace..."; \
	BASE_IMAGE=$$(grep '^ARG BASE_IMAGE=' Dockerfile | head -1 | cut -d= -f2); \
	PROFILE=$$(echo "$$BASE_IMAGE" | sed 's|.*/||; s|tidyverse|analysis|; s|verse|rocker/verse|; s|r-ver|minimal|'); \
	USERNAME=$$(grep '^ARG USERNAME=' Dockerfile | head -1 | cut -d= -f2); \
	USERNAME=$${USERNAME:-analyst}; \
	HOME_DIR="/home/$$USERNAME"; \
	echo "🐳 Starting R ($$PROFILE)..."; \
	echo ""; \
	mkdir -p $$HOME/.cache/R/renv 2>/dev/null || true; \
	$(CONTAINER_RUNTIME) run --rm -it \
		--platform $(PLATFORM) \
		-v $$(pwd):$$HOME_DIR/project \
		-v $$HOME/.cache/R/renv:/opt/renv/cache \
		-w $$HOME_DIR/project \
		-e KITTY_WINDOW_ID="$${KITTY_WINDOW_ID:-}" \
		-e ITERM_SESSION_ID="$${ITERM_SESSION_ID:-}" \
		-e TERM_PROGRAM="$${TERM_PROGRAM:-}" \
		-e GHOSTTY_RESOURCES_DIR="$${GHOSTTY_RESOURCES_DIR:-}" \
		-e WEZTERM_EXECUTABLE="$${WEZTERM_EXECUTABLE:-}" \
		$(PACKAGE_NAME) R || true; \
	echo ""; \
	echo "📋 Post-session validation..."; \
	if Rscript -e "quit(status = if(requireNamespace('zzrenvcheck', quietly=TRUE)) 0L else 1L)" 2>/dev/null; then \
		Rscript -e "zzrenvcheck::check_packages(auto_fix = TRUE, strict = TRUE)" \
		|| echo "⚠️  Package validation failed"; \
	else \
		echo "zzrenvcheck not installed on host."; \
		printf "Install now? [Y/n] "; \
		read _ans; \
		case "$$_ans" in \
			[Nn]*) echo "Skipping validation." ;; \
			*) Rscript -e "install.packages('remotes'); remotes::install_github('rgt47/zzrenvcheck', upgrade = 'never')" \
			   && Rscript -e "zzrenvcheck::check_packages(auto_fix = TRUE, strict = TRUE)" \
			   || echo "⚠️  Package validation failed" ;; \
		esac; \
	fi
	@if [ -f Dockerfile ] && $(CONTAINER_RUNTIME) image inspect $(PACKAGE_NAME) >/dev/null 2>&1; then \
		$(MAKE) test-container; \
	fi

docker-rstudio: ## RStudio Server at http://localhost:8787
	@echo "Starting RStudio Server on http://localhost:8787"
	@echo "Username: rstudio, Password: rstudio"
	@echo "Terminal available for code editing with vim"
	$(CONTAINER_RUNTIME) run --rm -it -p 8787:8787 -v $$(pwd):/home/rstudio/project $(PACKAGE_NAME) /init

# Alias for rstudio
rstudio: docker-rstudio ## Alias for docker-rstudio

##@ Native R (host renv, or Nix shell when a flake is present)

document: ## Generate roxygen2 documentation
	$(R_RUN) R --quiet -e "devtools::document()"

build: document ## R CMD build (documents first)
	$(R_RUN) R CMD build .

check: build ## R CMD check --as-cran (builds first)
	$(R_RUN) R CMD check --as-cran *.tar.gz

install: document ## devtools::install
	$(R_RUN) R --quiet -e "devtools::install()"

vignettes: document ## Build vignettes
	$(R_RUN) R --quiet -e "devtools::build_vignettes()"

test: ## Run the tinytest suite
	$(R_RUN) R --quiet -e "tinytest::run_test_dir('inst/tinytest')"

deps: ## Install package dependencies
	$(R_RUN) R --quiet -e "devtools::install_deps(dependencies = TRUE)"

##@ Docker build

# Pre-build audit: check renv.lock packages against PPM for missing system
# dependencies and packages that will compile from source.
# Requires python3 (standard on macOS/Linux). Pass SKIP_AUDIT=1 to bypass.
check-binaries: ## Audit renv.lock for missing system deps (pre-build)
	@if [ -f tools/check-binaries.py ]; then \
		python3 tools/check-binaries.py \
			--renv-lock renv.lock \
			--dockerfile Dockerfile; \
	else \
		echo "  tools/check-binaries.py not found, skipping audit"; \
	fi

docker-build: check-binaries ## Pre-build audit then build image from renv.lock
	zzcollab rebuild

docker-rebuild: ## Rebuild image without cache (force fresh build)
	zzcollab rebuild --no-cache

docker-build-log: ## Build image with detailed logs (for debugging)
	zzcollab rebuild --log

##@ Apptainer (HPC)

# Build the SIF from the project's local docker/podman image (build it first
# with `make docker-build`). On a cluster you would instead build from the
# digest-pinned team image: apptainer build $(SIF) docker://ACCOUNT/IMAGE@sha256:...
sif: ## Build env.sif from the project image (apptainer runtime)
	apptainer build $(SIF) docker-daemon://$(PACKAGE_NAME):latest
	@echo "✅ Built $(SIF). Use it: make CONTAINER_RUNTIME=apptainer r"

##@ Docker package operations

docker-document: ## roxygen2 documentation in the container
	$(DOCKER_RUN) R --quiet -e "devtools::document()"

docker-build-pkg: docker-document ## R CMD build in the container
	$(DOCKER_RUN) R CMD build .

docker-check: docker-build-pkg ## R CMD check in the container
	$(DOCKER_RUN) R CMD check --as-cran *.tar.gz

docker-test: ## Run tinytest in the container (uncached)
	$(DOCKER_RUN) R --quiet -e "tinytest::run_test_dir('inst/tinytest')"

# Cached tinytest target. Re-runs only when R sources, test files, or
# renv.lock change since the last passing run. Add .make/ to .gitignore.
TEST_SOURCES := $(wildcard R/*.R) $(wildcard inst/tinytest/*.R)

.make/tests-passed: $(TEST_SOURCES) renv.lock
	@mkdir -p .make
	$(CONTAINER_RUNTIME) run --rm \
		--platform $(PLATFORM) \
		-v "$(CURDIR)":$(CONTAINER_HOME)/project \
		-w $(CONTAINER_HOME)/project \
		-e ZZCOLLAB_CONTAINER=true \
		$(PACKAGE_NAME) \
		Rscript -e "if (dir.exists('inst/tinytest')) \
		  tinytest::run_test_dir('inst/tinytest') \
		  else message('No inst/tinytest -- skipping')"
	@touch $@

test-container: .make/tests-passed ## Run tinytest in the container (cached by make)

docker-vignettes: docker-document ## Build vignettes in the container
	$(DOCKER_RUN) R --quiet -e "devtools::build_vignettes()"

docker-render: ## Render REPORT (default analysis/report/report.Rmd)
	$(DOCKER_RUN) R --quiet -e "rmarkdown::render('$(REPORT)')"

# Run an arbitrary R script in the container (used by zzcollab::run_script)
docker-script:
	$(DOCKER_RUN) Rscript -e "source('$(SCRIPT)')"

docker-render-qmd: ## Render the Quarto report (analysis/report/index.qmd)
	$(DOCKER_RUN) quarto render analysis/report/index.qmd

# Back-compat aliases. docker-check-renv is invoked by zzcollab::run_script;
# both are identical to check-renv-no-fix and check-renv respectively, so they
# delegate rather than duplicate the docker run command.
docker-check-renv: check-renv-no-fix ## Validate deps in container (alias)

docker-check-renv-fix: check-renv

##@ Team images

# Multi-arch team publishing uses docker buildx + imagetools, which have no
# podman equivalent, so these targets are docker-specific regardless of
# CONTAINER_RUNTIME. Local build (zzcollab rebuild) and run honour the runtime.
docker-push-team: ## Multi-arch build and push; writes .team-image-digest
	@if [ -z "$(DOCKERHUB_ACCOUNT)" ]; then \
	  echo "DOCKERHUB_ACCOUNT is not set." >&2; \
	  echo "Run: make docker-push-team DOCKERHUB_ACCOUNT=youraccount" >&2; \
	  exit 1; \
	fi
	@echo "Building and pushing multi-arch image $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG)"
	@echo "Platforms: linux/amd64,linux/arm64 (native on both x86 and Apple Silicon)"
	docker buildx build \
		--platform linux/amd64,linux/arm64 \
		--tag $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG) \
		--push \
		.
	@echo "✅ Team image pushed: $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG)"
	@echo "   Recording image digest for reproducible pulls..."
	@docker buildx imagetools inspect \
		$(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG) \
		--format '{{json .Manifest.Digest}}' \
		| tr -d '"' > .team-image-digest
	@echo "   Digest written to .team-image-digest"
	@echo "   Commit .team-image-digest so team members pull the exact same image."

# Pull the team image pinned to the digest in .team-image-digest.
# This guarantees every team member runs bit-identical containers regardless
# of when they pull (the mutable tag may be overwritten; the digest never is).
docker-pull-team: ## Pull team image pinned by digest in .team-image-digest
	@if [ -z "$(DOCKERHUB_ACCOUNT)" ]; then \
	  echo "DOCKERHUB_ACCOUNT is not set." >&2; \
	  echo "Run: make docker-pull-team DOCKERHUB_ACCOUNT=youraccount" >&2; \
	  exit 1; \
	fi
	@if [ ! -f .team-image-digest ]; then \
	  echo "No .team-image-digest found. Run 'make docker-push-team' first." >&2; \
	  exit 1; \
	fi
	@DIGEST=$$(cat .team-image-digest); \
	REF="$(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME)@$$DIGEST"; \
	echo "Pulling team image pinned to digest: $$DIGEST"; \
	docker pull "$$REF"; \
	docker tag "$$REF" $(PACKAGE_NAME); \
	echo "✅ Team image loaded as $(PACKAGE_NAME)"

##@ Data integrity

# Write a sha256 manifest of the raw-data dir so any silent mutation is
# detectable. Commit data-manifest.sha256 and regenerate after data updates.
hash-data: ## Write sha256 manifest of the raw-data dir
	@if [ ! -d "$(RAW_DATA_DIR)" ]; then \
	  echo "No $(RAW_DATA_DIR) directory found; skipping hash-data."; \
	  exit 0; \
	fi
	find "$(RAW_DATA_DIR)" -type f | sort | xargs shasum -a 256 > data-manifest.sha256
	@echo "Data manifest written to data-manifest.sha256"
	@echo "Commit this file to make future mutations detectable."

verify-data: ## Verify data files against data-manifest.sha256
	@if [ ! -f data-manifest.sha256 ]; then \
	  echo "No data-manifest.sha256 found; run 'make hash-data' first." >&2; \
	  exit 1; \
	fi
	shasum -a 256 --check data-manifest.sha256
	@echo "All data files verified."

##@ Cleanup

clean: ## Remove build artifacts (*.tar.gz, *.Rcheck)
	rm -f *.tar.gz
	rm -rf *.Rcheck

docker-clean: ## Remove this project's image (system-wide: docker-prune-*)
	$(CONTAINER_RUNTIME) rmi $(PACKAGE_NAME) || true

docker-disk-usage: ## Show Docker disk usage
	@echo "Docker disk usage:"
	@docker system df

docker-prune-cache: ## Remove Docker build cache (system-wide)
	@echo "Removing Docker build cache..."
	docker builder prune -af
	@echo "✅ Build cache cleaned"
	@$(MAKE) docker-disk-usage

docker-prune-all: ## Deep clean all unused Docker resources (system-wide)
	@echo "WARNING: This will remove all unused Docker images, containers, and build cache"
	@echo "Press Ctrl+C to cancel, or Enter to continue..."
	@read dummy
	@echo "Removing all unused Docker resources..."
	docker system prune -af
	@echo "✅ Docker cleanup complete"
	@$(MAKE) docker-disk-usage

.PHONY: document build check install vignettes test deps snapshot reset-lock ensure-zzrenvcheck style lint check-renv check-renv-no-fix check-renv-no-strict check-binaries docker-build docker-rebuild docker-build-log docker-push-team docker-pull-team docker-document docker-build-pkg docker-check docker-test docker-vignettes docker-render docker-script docker-render-qmd docker-rstudio rstudio r test-container docker-check-renv docker-check-renv-fix hash-data verify-data clean docker-clean docker-disk-usage docker-prune-cache docker-prune-all help
