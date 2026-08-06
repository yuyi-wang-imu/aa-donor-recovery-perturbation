# Environment records

`requirements.txt` defines the supported installation ranges for the
non-Geneformer Python workflows. `python_packages.tsv` records the packages
observed on the Windows repository-audit host; it is provenance and must not be
treated as a lock file. Some audit-host versions are newer than the supported
workflow ranges.

The recorded GPU Geneformer environment is documented separately in
`geneformer_gpu_environment_20260804_v1.yml`. R package versions and the
installation helper are provided in `r_packages.tsv` and
`install_r_packages.R`. These records document the environments used or
checked; they do not bundle external packages, CUDA, model weights, or licensed
software.
