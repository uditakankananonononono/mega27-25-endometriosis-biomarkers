# Endometriosis - proposed separate project

Status: project scaffold and gate audit only, not a finished paper or a positive benchmark. The parent repository is shared core; this directory must hold a disease-specific protocol, accession and external-service evidence ledger, reproducible results and a 50-page substantive paper before its gates can be claimed.

Current disease-tagged manifest records: 101 = 14 GSE studies + 86 nested GSM samples + 1 other. These are record units, not independent datasets or patients. The shared 40-service and 49-page PDF do not transfer as automatic per-project passes. Benchmark/discovery endpoint is open.

Disease-specific documents and source logs are not yet split from shared core; use the shared source code and result filenames by disease as leads, then verify original record attribution before copying.

## GSE7846 tissue-label contradiction

The historically used [GSE7846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7846) discovery series contains five case and five control arrays. Its series title, study summary and all sample **titles** say cultured human endometrial endothelial cells (HEECs) from eutopic endometrium. But individual full GEO records have a systematic split: all five case **descriptions** say "Gene expression data from ectopic endometrium," while all five control descriptions say eutopic. Individual characteristics name endothelial cells from a patient with endometriosis or a normal control, without resolving the anatomical source. `sources/GSE7846_description_conflict.csv` preserves every fetched URL, response hash and competing fields. Thus case/control diagnosis is source-grounded, but the case tissue's eutopic-versus-ectopic identity is not, and this series must not support a lesion-origin versus eutopic design claim. It remains a previously analyzed record-level cohort with disclosed uncertainty; no ten GSMs are added to the used manifest here while the endometriosis evidence ledger is under review. A depositor clarification or patient-to-site record is needed before any anatomical inference. The affected endpoint stays open.
