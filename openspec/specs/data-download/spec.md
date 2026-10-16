# data-download Specification

## Purpose

Fetches SEA-AD files from the public `sea-ad-single-cell-profiling` S3 bucket into a configurable local data root, so every machine (laptop, office, home, cluster) resolves the same relative layout and downloads are idempotent.

## Requirements

### Requirement: Data root resolution
The system SHALL resolve the local data root in this precedence: an explicit command-line destination, then the `SEAD_DATA_ROOT` environment variable, then `data/` relative to the current working directory. Every command that reads or writes raw data SHALL use the same resolution and SHALL print the resolved root at start.

#### Scenario: Environment variable sets the root
- **WHEN** `SEAD_DATA_ROOT=/home1/tgonza/sead_sexdiff/data` is set and `download` runs without a destination flag
- **THEN** files are written under `/home1/tgonza/sead_sexdiff/data/<s3 key>` and the printed root matches that path

#### Scenario: Command-line destination overrides the environment
- **WHEN** both `SEAD_DATA_ROOT` and an explicit `--dest-root` are given
- **THEN** the explicit destination is used and the environment value is ignored

#### Scenario: Default root
- **WHEN** neither is set
- **THEN** the root is `./data`

### Requirement: MTG donor-object download per omic
The system SHALL download the MTG per-donor objects when an omic is selected explicitly: `rna` maps to prefix `MTG/RNAseq/donors_objects/`, `atac` maps to `MTG/ATACseq/donor_objects/`, and `both` fetches both. Files SHALL be stored under the data root preserving the full S3 key. Without an omic selection the existing prefix/suffix behaviour and its defaults SHALL remain unchanged.

#### Scenario: Download both omics
- **WHEN** `download --omic both` runs against the public bucket
- **THEN** 84 RNA `.h5ad` files and 86 ATAC `.h5ad` files exist under `<root>/MTG/RNAseq/donors_objects/` and `<root>/MTG/ATACseq/donor_objects/`

#### Scenario: Generic prefix still works
- **WHEN** `download --prefix MTG/ATACseq/ --file SEAAD_MTG_ATACseq_all-nuclei_metadata.2024-12-06.csv` runs
- **THEN** only that file (plus `Changelog.html` if present) is downloaded under the preserved key

#### Scenario: Existing default is untouched
- **WHEN** `download` runs with no `--omic` and no `--prefix`
- **THEN** it behaves as before, fetching CSV and HTML files under `PFC/RNAseq/`

### Requirement: Idempotent, size-checked downloads
The system SHALL skip any file already present locally whose size equals the remote object size, SHALL re-download files whose size differs, and SHALL report counts of downloaded and skipped files and the total bytes transferred.

#### Scenario: Re-running after a complete download
- **WHEN** `download --omic rna` runs a second time with all 84 files present at the correct size
- **THEN** no file is transferred and the report states 0 downloaded, 84 skipped

#### Scenario: Truncated file is repaired
- **WHEN** one local file is smaller than the remote object
- **THEN** that file is re-downloaded and the others are skipped

### Requirement: Download manifest
The system SHALL write a manifest (JSON) next to the downloaded prefix listing each file's S3 key, size, remote ETag and the download timestamp, so downstream steps can record which release files produced their outputs.

#### Scenario: Manifest after download
- **WHEN** a download of the ATAC donor objects completes
- **THEN** `<root>/MTG/ATACseq/donor_objects/download_manifest.json` lists 86 entries with key, size, etag and timestamp
