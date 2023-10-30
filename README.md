
This repository is dedicated to the analysis undergone by Max Bladen between August and October (2023) for Patrick Albertus and Donald Lynch from the University of Cincinnati. The repository is structued as follows:

- `Config`:
  - Contains a single `.R` file which handles the configuration of the `.Rmd` knitting process. Not relevant to actual analysis, purely to maintain consistency across reports
- `Markdowns`:
  - Series of `.Rmd` files, each corresponding to a single report file. All code can be found within these files. All analysis and figures should be reproducible by knitting these documents,
  - Note: some of these `.Rmd` files yield `.rds` files which downstream markdowns are dependent on.
- `RDS`:
  - Contains any `.rds` objects yielded and/or used by the `.Rmd` files in `Markdowns`. Subsequent analysis can be completed on these `.rds` objects
- `Reports`:
  - `.html` files yielded by knitting the `.Rmd` files in `Markdowns`. Best to download/clone this repository and view these `.html` files locally, as Github cannot render them properly
- `Scripts`:
  - Contains any `.R` files used as part of the analysis that are outside the `.Rmd`. Currently only contains `Functions.R` which is a sourced script used to store some custom functions which were repeatedly used.
