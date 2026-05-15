# CDSS - Chirp-DSSS-Polar digital modem build system.
#
# This Makefile is intentionally thin: FPM remains the authoritative Fortran
# build driver, while make provides reproducible project-level workflows for
# release builds, debug builds, tests, installation, and common CLI runs.
#
# External requirements:
#   - gfortran with Fortran 2018 support
#   - fpm
#   - FFTW3 development headers and library (fftw3.f03 and libfftw3)
#
# Override examples:
#   make build FPM=/path/to/fpm
#   make install PREFIX=/opt/cdss
#   make run ARGS='info'

# Fortran Package Manager executable. Kept configurable for CI and local tools.
FPM ?= fpm

# Installation prefix used by `make install`. The binary is copied to bin/cdss.
PREFIX ?= $(HOME)/.local

# Parallelism used by OpenMP-enabled release/test runs. `nproc` is Linux-specific,
# which matches the current project target environment.
NUM_CPUS := $(shell nproc)

# Release/debug/test builds all need FFTW's Fortran include file. The Makefile
# exports flags so `fpm test` can find fftw3.f03 even when FPM is run through make.
# The release-oriented flags favor executable throughput over strict IEEE debug
# behavior; use `make debug` while investigating numerical issues.
export FPM_FFLAGS  = -O3 -fopenmp -march=native -flto -I/usr/include
export FPM_LDFLAGS = -fopenmp -flto -lfftw3 -lm
export OMP_NUM_THREADS = $(NUM_CPUS)

# FPM stores binaries under compiler/profile-specific build directories. Select
# the newest executable so stale build directories do not shadow fresh builds.
FIND_BIN = find build -path "*/app/cdss" -executable -type f -printf '%T@ %p\n' | sort -nr | head -1 | cut -d' ' -f2-

.PHONY: all build debug test install info bench run clean help

# Default target: produce an optimized executable.
all: build

# Optimized build for normal command-line use.
build:
	@echo "[Build] Compiling CDSS in release mode ($(NUM_CPUS) OpenMP threads)..."
	@$(FPM) build --profile release

# Debug build for development. FPM controls the debug flags for this profile.
debug:
	@echo "[Build] Compiling CDSS in debug mode..."
	@$(FPM) build --profile debug

# Full regression suite. Tests intentionally print diagnostic measurements to
# make numerical or receiver regressions easier to inspect from CI logs.
test:
	@echo "[Test] Running CDSS test suite (debug profile)..."
	@$(FPM) test --profile debug

# Install the newest built executable to $(PREFIX)/bin/cdss.
install: build
	@mkdir -p $(PREFIX)/bin
	@bin=$$($(FIND_BIN)); \
	if [ -z "$$bin" ]; then echo "Error: binary not found. Run 'make build' first."; exit 1; fi; \
	cp "$$bin" $(PREFIX)/bin/cdss
	@echo "[Install] Installed cdss to $(PREFIX)/bin/cdss"

# Print modem defaults as JSON through the freshly built executable.
info: build
	@bin=$$($(FIND_BIN)); \
	if [ -z "$$bin" ]; then echo "Error: binary not found."; exit 1; fi; \
	$$bin info

# Internal throughput check using a deterministic diagnostic run.
bench: build
	@bin=$$($(FIND_BIN)); \
	if [ -z "$$bin" ]; then echo "Error: binary not found."; exit 1; fi; \
	echo "--- Single-thread ---"; \
	OMP_NUM_THREADS=1 $$bin simulate --bits 64 --snr-db 35 --seed 1; \
	echo "--- $(NUM_CPUS)-thread ---"; \
	OMP_NUM_THREADS=$(NUM_CPUS) $$bin simulate --bits 64 --snr-db 35 --seed 1

# Run arbitrary CLI arguments against the newest built executable.
run: build
	@bin=$$($(FIND_BIN)); \
	if [ -z "$$bin" ]; then echo "Error: binary not found."; exit 1; fi; \
	$$bin $(ARGS)

# Remove FPM build artifacts. FPM may ask for confirmation depending on version;
# piping "y" keeps the target non-interactive for CI and automation.
clean:
	@echo "[Clean] Removing FPM build artifacts..."
	@printf 'y\n' | $(FPM) clean

# Project-level command summary. Run `cdss help` after installation for CLI help.
help:
	@echo "CDSS project targets"
	@echo "--------------------"
	@echo "  make build      Compile optimized release binary"
	@echo "  make debug      Compile debug profile"
	@echo "  make test       Run full regression suite"
	@echo "  make install    Install cdss to PREFIX/bin (default: $(PREFIX)/bin)"
	@echo "  make info       Print modem geometry JSON"
	@echo "  make bench      Compare 1-thread and OpenMP diagnostic runs"
	@echo "  make run ARGS='info'"
	@echo "  make clean      Remove FPM build artifacts"
