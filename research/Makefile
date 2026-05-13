# Research and scientific-validation workflows for CDSS.
#
# This Makefile is intentionally kept under research/ so the repository root
# remains focused on building, testing, installing, and running the modem.
# Monte Carlo campaigns, Python dependencies, and publication-grade plots are
# research artifacts rather than normal end-user modem operations.

PYTHON ?= python3
VENV ?= .venv
PY := $(VENV)/bin/python
PIP := $(VENV)/bin/pip
REPO_ROOT := $(abspath ..)
RESULT_ROOT ?= mc_results
RUN_ID ?= simulate_mc_$(shell date +%Y%m%d_%H%M%S)
OUT_DIR ?= $(RESULT_ROOT)/$(RUN_ID)
LAST_RUN_FILE ?= .last_run
MPLCONFIGDIR ?= .matplotlib-cache

# Campaign controls. Override these on the command line, for example:
#   make quick WORKERS=8 OMP_THREADS=1
#   make publication PAYLOADS=32,64 SNRS=45,40,35 TRIALS_PER_POINT=500
BINARY ?=
PRESET ?= smoke
WORKERS ?= 4
OMP_THREADS ?= 1
SEED_BASE ?= 20260512
TIMEOUT_S ?= 240
PAYLOADS ?=
SNRS ?=
PROFILES ?=
TRIALS_PER_POINT ?=
INPUT ?=
PLOT_OUT ?=

.PHONY: help venv deps build-cdss run smoke quick publication plots analyze smoke-with-plots

help:
	@echo "CDSS research targets"
	@echo "---------------------"
	@echo "  make venv                 Create .venv and install plotting dependencies"
	@echo "  make smoke                Run a tiny Monte Carlo sanity campaign"
	@echo "  make quick                Run the development-sized campaign"
	@echo "  make publication          Run the full scientific campaign"
	@echo "  make plots                Generate CSV analysis tables and figures for .last_run"
	@echo "  make plots INPUT=path     Analyze a specific results.jsonl file"
	@echo "  make smoke-with-plots     Run smoke campaign and immediately plot it"
	@echo ""
	@echo "Useful overrides: WORKERS, OMP_THREADS, PAYLOADS, SNRS, PROFILES,"
	@echo "TRIALS_PER_POINT, OUT_DIR, INPUT, PLOT_OUT, BINARY."

$(VENV)/bin/python:
	@echo "[Research] Creating Python virtual environment at $(VENV)..."
	@$(PYTHON) -m venv $(VENV)

$(VENV)/.deps.stamp: requirements.txt | $(VENV)/bin/python
	@echo "[Research] Installing Python dependencies..."
	@$(PIP) install --upgrade pip
	@$(PIP) install -r requirements.txt
	@touch $@

venv: $(VENV)/.deps.stamp

deps: venv

build-cdss:
	@$(MAKE) -C $(REPO_ROOT) build

run: build-cdss venv
	@set -e; \
	bin="$(BINARY)"; \
	if [ -z "$$bin" ]; then \
		bin=$$(find $(REPO_ROOT)/build -path "*/app/cdss" -executable -type f -printf '%T@ %p\n' | sort -nr | head -1 | cut -d' ' -f2-); \
	fi; \
	if [ -z "$$bin" ]; then echo "Error: cdss binary not found. Run make -C .. build or set BINARY=..."; exit 2; fi; \
	extra_args=""; \
	if [ -n "$(PAYLOADS)" ]; then extra_args="$$extra_args --payloads $(PAYLOADS)"; fi; \
	if [ -n "$(SNRS)" ]; then extra_args="$$extra_args --snrs $(SNRS)"; fi; \
	if [ -n "$(PROFILES)" ]; then extra_args="$$extra_args --profiles $(PROFILES)"; fi; \
	if [ -n "$(TRIALS_PER_POINT)" ]; then extra_args="$$extra_args --trials-per-point $(TRIALS_PER_POINT)"; fi; \
	mkdir -p "$(OUT_DIR)"; \
	echo "[Research] Running Monte Carlo preset=$(PRESET) output=$(OUT_DIR)"; \
	$(PY) scripts/cdss_monte_carlo.py \
		--preset "$(PRESET)" \
		--binary "$$bin" \
		--out-dir "$(OUT_DIR)" \
		--workers "$(WORKERS)" \
		--omp-threads "$(OMP_THREADS)" \
		--seed-base "$(SEED_BASE)" \
		--timeout-s "$(TIMEOUT_S)" \
		$$extra_args; \
	printf '%s\n' "$(OUT_DIR)" > "$(LAST_RUN_FILE)"

smoke: PRESET=smoke
smoke: run

quick: PRESET=quick
quick: run

publication: PRESET=publication
publication: run

plots: venv
	@set -e; \
	input="$(INPUT)"; \
	if [ -z "$$input" ]; then \
		if [ ! -f "$(LAST_RUN_FILE)" ]; then echo "Error: no INPUT=... and no $(LAST_RUN_FILE). Run make smoke/quick/publication first."; exit 2; fi; \
		run_dir=$$(cat "$(LAST_RUN_FILE)"); \
		input="$$run_dir/results.jsonl"; \
	fi; \
	if [ ! -f "$$input" ]; then echo "Error: input JSONL not found: $$input"; exit 2; fi; \
	out_dir="$(PLOT_OUT)"; \
	if [ -z "$$out_dir" ]; then out_dir="$$(dirname "$$input")/analysis"; fi; \
	mkdir -p "$(MPLCONFIGDIR)"; \
	echo "[Research] Plotting $$input -> $$out_dir"; \
	MPLCONFIGDIR="$(MPLCONFIGDIR)" $(PY) scripts/plot_mc.py --input "$$input" --out-dir "$$out_dir"

analyze: plots

smoke-with-plots:
	@$(MAKE) smoke OUT_DIR="$(OUT_DIR)"
	@$(MAKE) plots INPUT="$(OUT_DIR)/results.jsonl" PLOT_OUT="$(OUT_DIR)/analysis"
