# Group adjacent peaks by non-overlapping windows

Group peaks into non-overlapping genomic windows of size `window_size`,
using peak coordinates from a BED-like metadata table.

## Usage

``` r
cacti_group_peak_window(
  window_size,
  file_pheno_meta,
  file_peak_group,
  file_peak_group_peaklevel
)
```

## Arguments

- window_size:

  Character like `"50kb"` (kb units) or numeric (bp).

- file_pheno_meta:

  Path to input phenotype meta table (with header). Columns are -

  - `phe_chr`: chromosome, with "chr" prefix (e.g. "chr1")

  - `phe_from`: peak start (integer, \< `phe_to`)

  - `phe_to`: peak end (integer)

  - `phe_id`: peak identifier (unique per row)

- file_peak_group:

  Output path for window-level table. Columns are - group, group_chr,
  group_from, group_to, pid_group, n_pid_group

- file_peak_group_peaklevel:

  Output path for peak-level table. Columns are - phe_id, phe_chr,
  phe_from, phe_to, group, n_pid_group

## Value

Invisibly returns the grouped window-level table; writes
`file_peak_group`.
