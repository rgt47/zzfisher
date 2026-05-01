# zzcollab Makefile v2.4.0
# Docker-first workflow for reproducible research

# Auto-detect from project (no manual configuration needed)
PACKAGE_NAME := $(shell basename $(CURDIR))
PROJECT_NAME := $(PACKAGE_NAME)
R_VERSION := $(shell grep 'R_VERSION=' Dockerfile 2>/dev/null | head -1 | sed 's/.*=//' || echo "4.4.0")

# Git-based versioning for reproducibility (use git SHA or date)
GIT_SHA := $(shell git rev-parse --short HEAD 2>/dev/null || echo "$(shell date +%Y%m%d)")
IMAGE_TAG = $(GIT_SHA)

# Help target (default)
help:
	@echo "Available targets:"
	@echo ""
	@echo "  Validation (NO HOST R REQUIRED!):"
	@echo "    check-renv             - Full validation: strict + auto-fix (recommended)"
	@echo "    check-renv-no-fix      - Validation only, no auto-install"
	@echo "    check-renv-no-strict   - Standard mode (skip tests/, vignettes/)"
	@echo "    check-system-deps      - Check for missing system dependencies in Dockerfile"
	@echo ""
	@echo "  Main workflow (RECOMMENDED):"
	@echo "    r                     - Start bash terminal (vim editing, all profiles)"
	@echo "    rstudio               - Start RStudio Server on http://localhost:8787"
	@echo ""
	@echo "  Native R - requires local R installation:"
	@echo "    document, build, check, install, vignettes, test, deps"
	@echo ""
	@echo "  Docker utilities:"
	@echo "    docker-build          - Build image from current renv.lock"
	@echo "    docker-rebuild        - Rebuild image without cache (force fresh build)"
	@echo "    docker-build-log      - Build with detailed logs (for debugging)"
	@echo "    docker-push-team, docker-document, docker-build-pkg, docker-check"
	@echo "    docker-test, docker-vignettes, docker-render, docker-render-qmd"
	@echo "    docker-check-renv"
	@echo ""
	@echo "  Cleanup:"
	@echo "    clean, docker-clean"
	@echo "    docker-prune-cache       - Remove Docker build cache"
	@echo "    docker-prune-all         - Deep clean (all unused Docker resources)"
	@echo "    docker-disk-usage        - Show Docker disk usage"

# Native R targets (require local R installation)
document:
	R --quiet -e "devtools::document()"

build:
	R CMD build .

check: document
	R CMD check --as-cran *.tar.gz

install: document
	R --quiet -e "devtools::install()"

vignettes: document
	R --quiet -e "devtools::build_vignettes()"

test:
	R --quiet -e "devtools::test()"

deps:
	R --quiet -e "devtools::install_deps(dependencies = TRUE)"

# Validate package dependencies (Pure shell, NO HOST R REQUIRED!)
# Checks that all packages used in code are in DESCRIPTION and renv.lock
# Full validation with strict mode, auto-fix, and verbose output (DEFAULT behavior)
# Scans all directories (root, R/, scripts/, analysis/, tests/, vignettes/, inst/)
# Auto-adds missing packages to DESCRIPTION and renv.lock
# Run this before `git commit` to catch issues locally (prevents CI failures)
check-renv:
	@zzcollab validate --fix --strict --verbose

# Validation only, no auto-fix (report issues without modifying files)
check-renv-no-fix:
	@zzcollab validate --no-fix --strict --verbose

# Standard mode validation (skip tests/, vignettes/, inst/ directories)
check-renv-no-strict:
	@zzcollab validate --fix --verbose

# Legacy: R-based validation (for CI/CD that has R pre-installed)
# This is the old approach, kept for backward compatibility
check-renv-ci:
	@zzcollab validate --fix --strict --verbose

# Check for missing system dependencies in Dockerfile
# Scans codebase for R packages and detects missing system libraries
# Suggests where to add missing deps in CUSTOM_SYSTEM_DEPS sections
check-system-deps:
	@zzcollab validate --system-deps

# Docker targets (work without local R)
# Docker-first workflow:
#   1. Work in containers (make r)
#   2. Install packages (renv::install("pkg"))
#   3. Exit container (auto-snapshot on exit)
#   4. Build new image (make docker-build)

docker-build:
	zzcollab build

docker-rebuild:
	zzcollab build --no-cache

docker-build-log:
	zzcollab build --log

docker-push-team:
	@echo "Tagging image as $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG)"
	docker tag $(PACKAGE_NAME) $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG)
	@echo "Pushing to Docker Hub..."
	docker push $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG)
	@echo "✅ Team image pushed: $(DOCKERHUB_ACCOUNT)/$(PROJECT_NAME):$(IMAGE_TAG)"
	@echo "   Team members should update .zzcollab_team_setup to reference this tag"

docker-document:
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R --quiet -e "devtools::document()"

docker-build-pkg:
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R CMD build .

docker-check: docker-document
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R CMD check --as-cran *.tar.gz

docker-test:
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R --quiet -e "devtools::test()"

docker-vignettes: docker-document
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R --quiet -e "devtools::build_vignettes()"

docker-render:
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R --quiet -e "rmarkdown::render('analysis/report/report.Rmd')"

docker-render-qmd:
	docker run --rm -v $$(pwd):/home/analyst/project -w /home/analyst/project $(PACKAGE_NAME) quarto render analysis/report/index.qmd

docker-check-renv:
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R --quiet -e "renv::status()"

docker-check-renv-fix:
	docker run --rm -v $$(pwd):/home/analyst/project $(PACKAGE_NAME) R --quiet -e "renv::snapshot()"

docker-rstudio:
	@echo "Starting RStudio Server on http://localhost:8787"
	@echo "Username: rstudio, Password: rstudio"
	@echo "Terminal available for code editing with vim"
	docker run --rm -it -p 8787:8787 -v $$(pwd):/home/rstudio/project $(PACKAGE_NAME) /init

# Terminal: Interactive bash for vim editing
r: check-renv
	@if [ ! -f Dockerfile ]; then \
		echo ""; \
		echo "❌ No Dockerfile found - workspace not initialized"; \
		echo ""; \
		echo "Run zzcollab to create a Docker environment:"; \
		echo ""; \
		echo "  zzcollab docker                            # default profile"; \
		echo "  zzcollab docker --profile analysis         # tidyverse"; \
		echo "  zzcollab docker --profile publishing       # with LaTeX"; \
		echo ""; \
		echo "See: zzcollab docker --help for all options"; \
		echo ""; \
		exit 1; \
	fi; \
	echo "🔍 Checking system dependencies..."; \
	zzcollab validate --system-deps 2>/dev/null || true; \
	BASE_IMAGE=$$(grep '^ARG BASE_IMAGE=' Dockerfile | head -1 | cut -d= -f2); \
	PROFILE=$$(echo "$$BASE_IMAGE" | sed 's|.*/||; s|tidyverse|analysis|; s|verse|publishing|; s|r-ver|minimal|'); \
	USERNAME=$$(grep '^ARG USERNAME=' Dockerfile | head -1 | cut -d= -f2); \
	USERNAME=$${USERNAME:-analyst}; \
	HOME_DIR="/home/$$USERNAME"; \
	echo "🐳 Starting R ($$PROFILE)..."; \
	echo ""; \
	mkdir -p $$HOME/.cache/R/renv 2>/dev/null || true; \
	docker run --rm -it \
		-v $$(pwd):$$HOME_DIR/project \
		-v $$HOME/.cache/R/renv:$$HOME_DIR/.cache/R/renv \
		-w $$HOME_DIR/project \
		-e KITTY_WINDOW_ID="$${KITTY_WINDOW_ID:-}" \
		-e ITERM_SESSION_ID="$${ITERM_SESSION_ID:-}" \
		-e TERM_PROGRAM="$${TERM_PROGRAM:-}" \
		-e GHOSTTY_RESOURCES_DIR="$${GHOSTTY_RESOURCES_DIR:-}" \
		-e WEZTERM_EXECUTABLE="$${WEZTERM_EXECUTABLE:-}" \
		$(PACKAGE_NAME) R; \
	echo ""; \
	echo "📋 Post-session validation..."; \
	zzcollab validate --fix --strict --verbose || echo "⚠️  Package validation failed"

# Alias for rstudio
rstudio: docker-rstudio

# Cleanup
clean:
	rm -f *.tar.gz
	rm -rf *.Rcheck

docker-clean:
	docker rmi $(PACKAGE_NAME) || true
	docker system prune -f

# Docker disk management
docker-disk-usage:
	@echo "Docker disk usage:"
	@docker system df

docker-prune-cache:
	@echo "Removing Docker build cache..."
	docker builder prune -af
	@echo "✅ Build cache cleaned"
	@make docker-disk-usage

docker-prune-all:
	@echo "WARNING: This will remove all unused Docker images, containers, and build cache"
	@echo "Press Ctrl+C to cancel, or Enter to continue..."
	@read dummy
	@echo "Removing all unused Docker resources..."
	docker system prune -af
	@echo "✅ Docker cleanup complete"
	@make docker-disk-usage

.PHONY: all document build check install vignettes test deps check-renv check-renv-no-fix check-renv-no-strict check-renv-ci docker-build docker-rebuild docker-build-log docker-push-team docker-document docker-build-pkg docker-check docker-test docker-vignettes docker-render docker-render-qmd docker-rstudio r docker-check-renv docker-check-renv-fix clean docker-clean docker-disk-usage docker-prune-cache docker-prune-all help
