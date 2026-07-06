# ==========================================
# zzcollab .Rprofile v0.1.0
# ==========================================
# Part 1: User Personal Settings (from ~/.Rprofile)
# Part 2: renv Activation + Reproducibility Options
# Part 3: Auto-Snapshot on Exit
# ==========================================

# ==========================================
# Part 1: User Personal Settings (always)
# ==========================================
q <- function(save="no", ...) quit(save=save, ...)

# Package installation behavior (non-interactive)
options(
  install.packages.check.source = "no",
  install.packages.compile.from.source = "never",
  Ncpus = parallel::detectCores()
)

# ==========================================
# RNG discipline (R-7)
# ==========================================
# Pin the RNG algorithm and normal-variate method explicitly so that
# stochastic analyses (bootstrap, MCMC, cross-validation, simulation)
# are reproducible across R versions. R 3.6.0 changed the default
# sample.kind, which silently breaks previously reproducible seeds.
# Set a project-level seed with set.seed() in each analysis script.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ==========================================
# Part 2: Container Detection
# ==========================================
# Set ZZCOLLAB_CONTAINER=true in Dockerfile to enable renv
in_container <- Sys.getenv("ZZCOLLAB_CONTAINER") == "true"

# Install mode governs whether renv manages packages. The Dockerfile records it
# (ZZCOLLAB_INSTALL_MODE) so this profile self-adapts when renv is toggled: in
# DESCRIPTION-install mode the renv workflow (auto-init, restore, snapshot) is
# skipped entirely, so removing renv does not get undone by container startup.
# Default "renv" so images predating this variable keep the full workflow.
install_mode <- Sys.getenv("ZZCOLLAB_INSTALL_MODE", "renv")
renv_enabled <- in_container && install_mode != "description"

# Set repos based on environment
if (in_container) {
  # Use Posit Package Manager for pre-compiled binaries in container. Derive the
  # mirror from RENV_CONFIG_REPOS_OVERRIDE, which the Dockerfile sets to the
  # dated PPM URL, so the Ubuntu codename and PPM snapshot live in exactly one
  # place (the Dockerfile) and cannot drift or be left unsubstituted. Set both
  # repos AND renv.repos.cran (renv uses the latter as its default).
  ppm_repo <- Sys.getenv("RENV_CONFIG_REPOS_OVERRIDE",
                         "https://packagemanager.posit.co/cran/__linux__/noble/latest")
  options(repos = c(CRAN = ppm_repo), renv.repos.cran = ppm_repo)
} else {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!renv_enabled) {
  # ==========================================
  # renv not active (host R, or container in DESCRIPTION-install mode)
  # ==========================================
  if (!in_container) {
    message("ℹ️ Host R session (renv skipped - use container for reproducibility)")
  } else {
    message("📦 Container R session (DESCRIPTION-install mode; renv not in use)")
  }

} else {
  # ==========================================
  # Container R: Full renv workflow
  # ==========================================

  message("🐳 Container R session (", Sys.getenv("HOSTNAME", "zzcollab"), ")")

  # CI detection (GitHub Actions sets CI=true)
  in_ci <- nzchar(Sys.getenv("CI"))

  # renv Cache Path Configuration
  # If RENV_PATHS_CACHE already set (e.g., via docker -e), use it
  # Otherwise use ~/.cache/R/renv (shared across projects)
  if (Sys.getenv("RENV_PATHS_CACHE") == "") {
    Sys.setenv(RENV_PATHS_CACHE = file.path(Sys.getenv("HOME"), ".cache/R/renv"))
  }

  # Activate renv (set project-local library paths)
  if (file.exists("renv/activate.R")) {
    source("renv/activate.R")
  } else if (in_container && nzchar(Sys.getenv("RENV_PATHS_LIBRARY"))) {
    # Image-library-authoritative mode: renv/ is not bind-mounted, so there
    # is no activate.R to source at runtime. The packages were baked into
    # RENV_PATHS_LIBRARY at build time; put that library on the path directly
    # so the project's declared packages resolve. renv's library layout is
    # <RENV_PATHS_LIBRARY>/<platform>/R-<major.minor>/<arch>.
    renv_lib_root <- Sys.getenv("RENV_PATHS_LIBRARY")
    baked_lib <- Sys.glob(file.path(renv_lib_root, "*", "R-*", "*"))
    baked_lib <- baked_lib[dir.exists(baked_lib)]
    if (length(baked_lib) > 0) {
      .libPaths(c(baked_lib[[1]], .libPaths()))
    } else {
      warning("⚠️  RENV_PATHS_LIBRARY set but no baked library found under ",
              renv_lib_root, call. = FALSE)
    }
  }

  # renv consent (skips first-time prompts)
  options(
    renv.consent = TRUE,
    renv.config.install.prompt = FALSE,
    renv.config.auto.snapshot = FALSE
  )

  # Helper function for initializing renv without prompts
  renv_init_quiet <- function() {
    renv::init(
      bare = TRUE,
      settings = list(snapshot.type = "implicit"),
      force = TRUE,
      restart = FALSE,
      load = FALSE
    )

    message("✅ renv initialized")
    message("   Install packages with: install.packages('package')")
  }

  # ==========================================
  # Auto-Initialize renv (New Projects)
  # ==========================================
  if (!file.exists("renv.lock")) {
    auto_init <- Sys.getenv("ZZCOLLAB_AUTO_INIT", "true")
    is_project <- file.exists("DESCRIPTION") || getwd() == "/home/analyst/project"

    if (tolower(auto_init) %in% c("true", "t", "1") && is_project) {
      message("\n🔧 ZZCOLLAB: Auto-initializing renv for new project...")
      tryCatch({
        renv_init_quiet()
      }, error = function(e) {
        warning("⚠️  Auto-init failed: ", conditionMessage(e),
                "\n   Run manually: renv_init_quiet()", call. = FALSE)
      })
    }
  } else {
    # ==========================================
    # Recover renv infrastructure if missing
    # ==========================================
    # This handles: renv.lock exists but renv/ doesn't (e.g., git clone on host).
    # Skip in the container: the baked library at RENV_PATHS_LIBRARY is the
    # source of truth; renv/ is not bind-mounted and recovery would attempt to
    # write to a read-only path.
    if (!file.exists("renv/activate.R") && !in_container) {
      message("\n🔧 ZZCOLLAB: renv.lock found but renv/ missing - recovering...")
      tryCatch({
        renv_init_quiet()
        if (file.exists("renv/activate.R")) {
          source("renv/activate.R")
        }
      }, error = function(e) {
        warning("⚠️  renv recovery failed: ", conditionMessage(e), call. = FALSE)
      })
    }

    # ==========================================
    # Auto-Restore Missing Packages
    # ==========================================
    auto_restore <- Sys.getenv("ZZCOLLAB_AUTO_RESTORE", "true")

    if (tolower(auto_restore) %in% c("true", "t", "1") && !in_ci) {
      in_lsp <- !interactive() || nzchar(Sys.getenv("NVIM_LISTEN_ADDRESS")) ||
                nzchar(Sys.getenv("RSTUDIO"))

      tryCatch({
        if (in_lsp) {
          invisible(suppressMessages(suppressWarnings({
            sink("/dev/null", type = "output")
            sink("/dev/null", type = "message")
            on.exit({
              sink(type = "message")
              sink(type = "output")
            }, add = TRUE)
            renv::restore(prompt = FALSE)
            sink(type = "message")
            sink(type = "output")
          })))
        } else {
          renv::restore(prompt = FALSE)
        }
      }, error = function(e) {
        if (!in_lsp && !grepl("already synchronized|consistent state", conditionMessage(e))) {
          warning("⚠️  Auto-restore failed: ", conditionMessage(e), call. = FALSE)
        }
        invisible(NULL)
      })
    }
  }

  # ==========================================
  # Auto-Snapshot on R Exit (Container only)
  # ==========================================
  .Last <- function() {
    auto_snapshot <- Sys.getenv("ZZCOLLAB_AUTO_SNAPSHOT", "true")

    if (tolower(auto_snapshot) %in% c("true", "t", "1") && !in_ci) {
      if (file.exists("renv.lock") && file.exists("renv/activate.R")) {
        message("\n📸 Auto-snapshot: Updating renv.lock...")

        snapshot_result <- tryCatch({
          renv::snapshot(prompt = FALSE)
          TRUE
        }, error = function(e) {
          warning("Auto-snapshot failed: ", conditionMessage(e), call. = FALSE)
          FALSE
        })

        if (snapshot_result) {
          message("✅ renv.lock updated successfully")
          message("   Commit changes: git add renv.lock && git commit -m 'Update packages'")
        }
      }
    }

    if (exists(".Last.user", mode = "function", envir = .GlobalEnv)) {
      tryCatch(
        .Last.user(),
        error = function(e) warning("User .Last failed: ", conditionMessage(e))
      )
    }
  }

  # Re-apply Posit PM repos AFTER renv::load() (which overrides from lockfile)
  options(repos = c(CRAN = ppm_repo))
}

# ==========================================
# Part 3: Reproducibility Options (always)
# ==========================================
# These ensure consistent behavior on both host and container
options(
  stringsAsFactors = FALSE,
  contrasts = c("contr.treatment", "contr.poly"),
  na.action = "na.omit",
  digits = 7,
  OutDec = "."
)

# ==========================================
# Part 4: Personal Customizations (always)
# ==========================================
if (file.exists(".Rprofile.local")) {
  tryCatch(
    source(".Rprofile.local"),
    error = function(e) {
      warning(".Rprofile.local failed to load: ", conditionMessage(e))
    }
  )
}
