# HiFi heterozygosity — Nextflow pipeline

DSL2 Nextflow port of `/data/het_L681/run_het_L681.sh`. Same six steps,
same `hifi-het:latest` Docker image, but with native parallel fan-out
across samples, `-resume` after failures, and an auto-generated DAG /
timeline / resource report.

## TL;DR

```bash
cd /data/het_L681/nextflow
./run.sh                       # runs all samples in samplesheet.csv
```

That's it — `run.sh` installs Nextflow + OpenJDK 17 on first run, then
executes `main.nf` with the default profile.

## Prerequisites

| Requirement | How to satisfy |
|---|---|
| Docker + `hifi-het:latest` image | `cd /data/het_L681 && ./00_setup.sh` |
| Java 17 (for Nextflow) | Auto-installed by `run.sh` if missing |
| Nextflow ≥ 23.x | Auto-installed by `run.sh` into `./nextflow` if missing |
| Input files exist on disk | Listed in [`samplesheet.csv`](samplesheet.csv) |

The pipeline does **not** rebuild the Docker image — `00_setup.sh` must
have run at least once on this host first.

## Files

| File | Purpose |
|---|---|
| [`main.nf`](main.nf) | Workflow + 5 processes: `INDEX_REF`, `MAP_HIFI`, `CALL_VARIANTS`, `FILTER_HET`, `COUNT_HET` |
| [`nextflow.config`](nextflow.config) | Docker on, default resources, `standard` + `parallel` profiles, timeline/report/DAG/trace enabled |
| [`samplesheet.csv`](samplesheet.csv) | One row per sample: `sample,assembly,reads` |
| [`run.sh`](run.sh) | Bootstraps Java + Nextflow, runs `main.nf`, forwards extra flags |

## Workflow

```
samplesheet.csv
       │  (one tuple per row)
       ▼
 ┌─────────────┐    ┌──────────┐    ┌───────────────┐    ┌────────────┐    ┌───────────┐
 │ INDEX_REF   │──▶ │ MAP_HIFI │──▶ │ CALL_VARIANTS │──▶ │ FILTER_HET │──▶ │ COUNT_HET │
 │ samtools    │    │ minimap2 │    │ bcftools      │    │ bcftools   │    │ bedtools  │
 │ faidx       │    │ map-hifi │    │ mpileup|call  │    │ view -i    │    │ intersect │
 └─────────────┘    └──────────┘    └───────────────┘    └────────────┘    └───────────┘
```

Each process runs inside the `hifi-het:latest` container with the
invoking user's UID/GID (so outputs aren't root-owned).

## Common invocations

```bash
# Default — local executor, one sample at a time
./run.sh

# Fan out across samples on a big instance (e.g. c7g.12xlarge)
./run.sh -profile parallel

# Resume after failure or after tuning a filter; cached steps are skipped
./run.sh -resume

# Override variant filters or window size on the CLI
./run.sh --window 50000 --min_qual 30 --min_dp 3 --max_dp 150

# Point at a different samplesheet
./run.sh --samplesheet /path/to/other.csv

# Custom output directory (default: ../results, i.e. /data/het_L681/results)
./run.sh --outdir /data/het_all/results
```

All flags after `./run.sh` are forwarded verbatim to `nextflow run main.nf`.

## Parameters

Defaults are set at the top of [`main.nf`](main.nf) and overridable
with `--<name> <value>`:

| Param | Default | Meaning |
|---|---|---|
| `--samplesheet` | `./samplesheet.csv` | CSV with columns `sample,assembly,reads` |
| `--outdir` | `../results` | Per-sample output goes to `<outdir>/<sample>/` |
| `--window` | `100000` | Window size (bp) for per-window het counts |
| `--min_qual` | `20` | bcftools `QUAL` floor in `FILTER_HET` |
| `--min_dp` | `5` | Per-site `FORMAT/DP` floor (drop if coverage is low: try 3) |
| `--max_dp` | `200` | Per-site `FORMAT/DP` ceiling — set to ≈ `3 × median_coverage` to suppress collapsed-repeat het inflation |

To pick a real `MAX_DP` from the data, run once with the default, then
inspect the `INFO/DP` distribution in `results/<sample>/<sample>.vcf.gz`
(see the `bcftools query` snippet in the parent README), then re-run
with `--max_dp <value> -resume` — only `FILTER_HET` and `COUNT_HET`
re-execute.

## Samplesheet

[`samplesheet.csv`](samplesheet.csv) — header + one row per sample:

```
sample,assembly,reads
L681,/data/genomes/Dlaeve_L681_Schlossteich2_Euphallic.p_ctg.fasta,/data/genomes/L681.fastq.gz
L451,/data/genomes/Dlaeve_L451_Glasgow_Aphallic.p_ctg.fasta,/data/genomes/L451.fastq.gz
```

Paths must be readable by the container — anything under `/data` works
because the host mounts it (Nextflow auto-mounts the paths it sees in
input channels).

## Outputs

```
results/
├── <sample>/
│   ├── <sample>.sorted.bam            (+ .bai, .flagstat.txt)
│   ├── <sample>.vcf.gz                (+ .csi)        — all variant sites
│   ├── <sample>.het.vcf.gz            (+ .csi)        — het SNPs only
│   ├── windows_100000.bed
│   ├── het_per_window.annotated.tsv   — per-window het counts + rate
│   └── het_per_contig.tsv             — per-contig totals, sorted by length
└── pipeline_info/
    ├── timeline.html   — Gantt of process executions
    ├── report.html     — CPU / RAM / IO per process
    ├── dag.html        — workflow graph
    └── trace.tsv       — machine-readable per-task metrics
```

Output schema for `het_per_window.annotated.tsv` and `het_per_contig.tsv`
is identical to the bash pipeline — see the parent README's "Output
schema" section.

## Profiles

Defined in [`nextflow.config`](nextflow.config):

| Profile | When | Behaviour |
|---|---|---|
| `standard` (default) | Single-sample run on a modest host | `local` executor, default resource hints |
| `parallel` | Big instance, many samples | `local` executor with executor-level `cpus = all`, `memory = 90 GB`; samples run concurrently as resources allow |

Process resource hints in `main.nf`:
- `MAP_HIFI`: up to 16 cpus, 24 GB RAM
- `CALL_VARIANTS`: up to 8 cpus, 8 GB RAM
- Others: config defaults (4 cpus, 8 GB RAM)

## Resume semantics

Nextflow caches each task by hash of (script + inputs + params). After
a failure or after tuning a parameter, `-resume` will re-execute only
the downstream tasks whose hash changed. Concretely:

- Add a new row to `samplesheet.csv` + `-resume` → only the new sample
  runs end-to-end; existing samples are skipped.
- Change `--max_dp` + `-resume` → `INDEX_REF`, `MAP_HIFI`, `CALL_VARIANTS`
  reuse cache; `FILTER_HET` + `COUNT_HET` re-run.
- Change `--window` + `-resume` → only `COUNT_HET` re-runs.

The cache lives in `./work/` — delete it to force a full rerun.

## Troubleshooting

- **`docker: permission denied`** — your user isn't in the `docker`
  group yet. Either re-login after `00_setup.sh`, run `newgrp docker`,
  or prefix commands with `sudo`.
- **Outputs owned by root** — shouldn't happen (the config sets
  `-u $(id -u):$(id -g)`). If it does, check `nextflow.config` wasn't
  overridden.
- **`hifi-het:latest` not found** — run `/data/het_L681/00_setup.sh`
  to build the image.
- **OOM in `MAP_HIFI`** — minimap2 + sort needs ~24 GB for a 1.8 Gb
  assembly. Lower `cpus` (less sort parallelism = lower peak RAM) or
  use a larger instance.
- **Slow `CALL_VARIANTS`** — bcftools mpileup is single-threaded per
  region; the `--threads` flag only parallelises compression. Expected.

## When to use this vs the bash pipeline

| | Bash (`../run_het_L681.sh`) | Nextflow (this folder) |
|---|---|---|
| One sample, one-off | ✓ simpler | works fine |
| 9 samples | manual loop | ✓ native fan-out (`-profile parallel`) |
| Resume after re-tuning filters | re-run from scratch | ✓ `-resume` |
| Reproducibility artefacts | scripts + Docker | scripts + Docker + DAG/timeline/report |
