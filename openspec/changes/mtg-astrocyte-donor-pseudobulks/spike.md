# Spike: MTG release contents (2026-10-16)

Files (downloaded with `download --prefix ... --file ...`, local size equals remote `ContentLength`):

| omic | key | size (B) | ETag |
|---|---|---|---|
| RNA | `MTG/RNAseq/donors_objects/H19.33.004_SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad` | 697,489,275 | `c6410739e76860b919ed03f0ac88f4d5-84` |
| ATAC | `MTG/ATACseq/donor_objects/H19.33.004_SEAAD_MTG_ATACseq_final-nuclei.2024-12-06.h5ad` | 79,165,509 | `62a4e8518384a4aa99c648537107b51b` |

Neither prefix holds a `Changelog.html` or any non-`.h5ad` object. Every file name matches
`<donor>_SEAAD_MTG_<RNAseq|ATACseq>_final-nuclei.<date>.h5ad`; one date per release
(RNA `2026-06-22`, ATAC `2024-12-06`).

Donor lists from the bucket: 84 RNA, 86 ATAC, every RNA donor has an ATAC object.
ATAC-only donors: `H19.30.002`, `H19.30.004`.

## obs columns

Both releases carry `Subclass`, `Supertype`, `Donor ID`, `Sex`,
`Overall AD neuropathological Change` (plus `Class`, `Genes detected`, `Used in analysis`, …).
ATAC additionally has `Supertype (non-expanded)` and `* confidence` columns. **D7 fallback is
not needed**: ATAC `obs` carries Sex and ADNC. Spike donor: `Female`, `Not AD` in both.

## Astrocyte labels

- `Subclass`: `Astrocyte` in **both** releases (the `Immune` vs `Microglia-PVM` difference does
  not affect astrocytes).
- `Supertype` among astrocytes, spike donor:
  - RNA: `Astro_1` 273, `Astro_2` 1465, `Astro_3` 400, `Astro_4` 14, `Astro_5` 90,
    `Astro_6-SEAAD` 24 (total 2266).
  - ATAC: `Astro_1` 4, `Astro_2` 20, `Astro_3` 1, `Astro_5` 9 (total 34).
- Mapping: identity for `Astro_1`…`Astro_5` and `Astro_6-SEAAD` in both releases (same
  taxonomy names). Any other label among astrocyte nuclei fails the run.

## Count matrices

- RNA: `X` is normalized float64 (not counts). Layer `UMIs` exists, CSR float64, all integer
  values (min nonzero 1, max 9337). Shape 27,431 × 36,601. → source `layers/UMIs`.
- ATAC: no layers; `X` is CSR float32 with integer values (min nonzero 1, max 378).
  Shape 1,540 × 218,882 (peak count matches the spec). → source `X`.
- Both are cast to `int64` after the integer check.

## var

- RNA: index = gene symbol; columns `gene_ids` (Ensembl), `feature_types`, `genome`.
- ATAC: index = peak `chrN:start-end`; no columns. Peaks are not in genomic order; order is
  kept as stored.

## QC notes

- RNA `Genes detected` equals the per-nucleus nnz of `UMIs` exactly. All 2266 astrocytes of
  the spike donor pass `min_genes=200` (min 1506 genes), so the RNA QC is essentially a no-op
  on `final-nuclei` objects; it is kept as specified.
- `Used in analysis` is `True` for every nucleus of the spike donor; not used as a filter.

## Reading strategy

anndata 0.12.7: `h5py` + `anndata.io.read_elem(f["obs"])` and
`anndata.io.sparse_dataset(f["layers/UMIs"])[rows]` reads only the astrocyte rows, so the
normalized RNA `X` (~1.7 GB in memory for this donor) is never loaded.

## Open questions answered

- Astrocyte `Subclass`/`Supertype` strings per release → above.
- ATAC matrix location → `X`, float32 integer-valued.
- ATAC `obs` carries Sex and ADNC → yes.
- `var` columns per release → above; copied through as-is.
