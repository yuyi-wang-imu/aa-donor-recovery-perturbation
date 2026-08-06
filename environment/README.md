# Environment records

`requirements.txt` defines the supported installation ranges for the
non-Geneformer Python workflows. `python_packages.tsv` records the packages
observed on the Windows repository-audit host; it is provenance and must not be
treated as a lock file. Some audit-host versions are newer than the supported
workflow ranges.

`publication_replay_python_20260806.txt` pins the Windows Python renderer used
for the clean publication-replay test. The replay verifier intentionally uses
documented pixel tolerances for Figures 3, 7 and S8. The frozen Figure 7 PNG can
also be reproduced byte-for-byte with the separately recorded
`figure7_exact_renderer_python_20260806.txt` stack; this alternate stack is not
required for scientific-value verification.
`publication_replay_r_packages_20260806.tsv` records the R 4.5.2 rendering
stack used in the successful clean replay.

The recorded GPU Geneformer environment is documented separately in
`geneformer_gpu_environment_20260804_v1.yml`. R package versions and the
installation helper are provided in `r_packages.tsv` and
`install_r_packages.R`. These records document the environments used or
checked; they do not bundle external packages, CUDA, model weights, or licensed
software.
