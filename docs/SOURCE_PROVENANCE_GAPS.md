# Source provenance and rights documentation gaps (endometriosis child, 2026-10-08 audit record)

This file records what this repository's own documents say about source hashes and terms. It is a documentation record, not a license verdict, not an integrity certificate, and not a clearance of any recorded rights hold. A matching claim below means two recorded hash strings are equal at a named field in pinned repository documents. It does not mean any assay bytes were re-downloaded or verified. No rights record was found for these sources, so reuse and redistribution are NOT cleared by this file.

Audited child commit: `28cc0cd504b62c9868fa2267263d515c031af003`  
Compared shared-parent repo: `mega27-25-biomarkers-underserved-diseases` at `36ff87995f4fe0e1d08f9e2b2cd885d341b88937`

Rights: No source-specific rights check found in inspected disease-child documentation. Hash matches are not legal clearance or byte verification.

## 1. Result-file hash records in this child (4)

Hash values are shown as 12-hex display prefixes only; the full value is in the cited file and field. Classes are kept separate on purpose: an expression/assay source hash, a GEO series-matrix hash (metadata that may include expression), a reference annotation, and a published artifact are different kinds of record.

- `results/endo_p40_result.json` `/sha256/GSE153739_GT_SO_5238_Transcript_Expression_Matrix.txt.gz` hash `627634eb5df0...` (GSE153739); class: non_source_payload_hash; equal to shared-parent field
- `results/endo_p40_result.json` `/sha256/gencode.v48.annotation.gtf.gz` hash `37a298c2b57e...` (GSE153739); class: reference_annotation_not_assay; equal to shared-parent field
- `results/endo_p17_result.json` `/source_sha256` hash `46b0404d797b...` (GSE212787); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/endo_p34_result.json` `/sha256` hash `e4d5dcad9c17...` (GSE313775); class: deposited_expression_source_file_record; equal to shared-parent field

## 4. Metadata crosswalk tables (not assay bytes)

These per-sample tables carry hash columns describing GEO source-response metadata. They are not blanket expression-matrix integrity.

- `projects/endometriosis/sources/GSE47360_used_sample_crosswalk.csv`: columns ['sha256'], 9 records; file identical to a shared-parent file: False; scope: GEO sample/series source-response metadata; not blanket expression matrix integrity

## Coverage boundary

- Source accessions listed in the child manifest: 14.
- Of those, 11 have no mapped result-file payload hash in this audit: GSE120103, GSE23339, GSE248593, GSE25628, GSE35287, GSE47360, GSE51981, GSE58178, GSE6364, GSE73622, GSE7846.
- 0 hash fields inherited from other diseases' records are not counted toward this child.
- Records absent from child manifest can still have result-file hashes. Neither presence nor absence proves assay acquisition by this scout. Other-disease copied hashes never count toward this child.
