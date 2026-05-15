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

# Campaign controls for the single release validation flow.
BINARY ?=
PRESET ?= release
WORKERS ?= 4
OMP_THREADS ?= 1
SEED_BASE ?= 6703739
TIMEOUT_S ?= 240
PAYLOADS ?=
SNRS ?=
PROFILES ?=
TRIALS_PER_POINT ?=
EXTRA_ARGS ?=
MIXED_SWEEP_ARGS ?=
INPUT ?=
PLOT_OUT ?=

# PDF report controls.
#   SKIP_BUILD=1        — skip fpm build (use when binary is already built)
#   REPORT_VERSION=v1.2 — embedded in filename and title page
#   REPORT_OUT=foo.pdf  — override default output path
#   REPORT_TITLE="..."  — override report title
SKIP_BUILD     ?= 0
REPORT_VERSION ?= dev
REPORT_TITLE   ?= CDSS Modem Validation Report
REPORT_OUT     ?=

.PHONY: help venv deps build-cdss run release report plots

help:
	@echo "CDSS research targets"
	@echo "---------------------"
	@echo "  make venv                 Create .venv and install plotting dependencies"
	@echo "  make release                 Run release validation and generate PDF"
	@echo "  make report                  Generate PDF report from last run (or INPUT=)"
	@echo "  make plots                   Generate analysis plots from last run (or INPUT=)"
	@echo ""
	@echo "Useful overrides: WORKERS, OMP_THREADS, OUT_DIR, INPUT, BINARY,"
	@echo "SKIP_BUILD, REPORT_VERSION, REPORT_OUT, REPORT_TITLE, PLOT_OUT."

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

run: venv
	@set -e; \
	if [ "$(PRESET)" != "release" ]; then \
		echo "Error: only PRESET=release is supported in this workflow."; exit 2; \
	fi; \
	if [ "$(SKIP_BUILD)" != "1" ]; then $(MAKE) -C $(REPO_ROOT) build; fi; \
	bin="$(BINARY)"; \
	if [ -z "$$bin" ]; then \
		bin=$$(find $(REPO_ROOT)/build -path "*/app/cdss" -executable -type f -printf '%T@ %p\n' | sort -nr | head -1 | cut -d' ' -f2-); \
	fi; \
	if [ -z "$$bin" ]; then echo "Error: cdss binary not found. Run make -C .. build or set BINARY=..."; exit 2; fi; \
	extra_args=""; \
	if [ -n "$(EXTRA_ARGS)" ]; then extra_args="$$extra_args $(EXTRA_ARGS)"; fi; \
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

release: PRESET=release
release: run report

report: venv
	@set -e; \
	input="$(INPUT)"; \
	if [ -z "$$input" ]; then \
		if [ ! -f "$(LAST_RUN_FILE)" ]; then \
			echo "Error: no INPUT=... and no $(LAST_RUN_FILE). Run a campaign first."; exit 2; \
		fi; \
		run_dir=$$(cat "$(LAST_RUN_FILE)"); \
		input="$$run_dir/results.jsonl"; \
	fi; \
	if [ ! -f "$$input" ]; then echo "Error: results.jsonl not found: $$input"; exit 2; fi; \
	out="$(REPORT_OUT)"; \
	if [ -z "$$out" ]; then \
		out_dir=$$(dirname "$$input"); \
		out="$$out_dir/cdss-$(REPORT_VERSION)-release-report.pdf"; \
	fi; \
	mkdir -p "$(MPLCONFIGDIR)"; \
	echo "[Research] Generating PDF report -> $$out"; \
	MPLCONFIGDIR="$(MPLCONFIGDIR)" MPLBACKEND=Agg \
	$(PY) scripts/generate_report.py \
		--input   "$$input" \
		--output  "$$out" \
		--version "$(REPORT_VERSION)" \
		--title   "$(REPORT_TITLE)"

plots: venv
	@set -e; \
	input="$(INPUT)"; \
	if [ -z "$$input" ]; then \
		if [ ! -f "$(LAST_RUN_FILE)" ]; then \
			echo "Error: no INPUT=... and no $(LAST_RUN_FILE). Run a campaign first."; exit 2; \
		fi; \
		run_dir=$$(cat "$(LAST_RUN_FILE)"); \
		input="$$run_dir/results.jsonl"; \
	fi; \
	if [ ! -f "$$input" ]; then echo "Error: input JSONL not found: $$input"; exit 2; fi; \
	out_dir="$(PLOT_OUT)"; \
	if [ -z "$$out_dir" ]; then out_dir="$$(dirname "$$input")/analysis"; fi; \
	mkdir -p "$(MPLCONFIGDIR)"; \
	echo "[Research] Plotting $$input -> $$out_dir"; \
	MPLCONFIGDIR="$(MPLCONFIGDIR)" MPLBACKEND=Agg \
	$(PY) scripts/plot_mc.py --input "$$input" --out-dir "$$out_dir"
