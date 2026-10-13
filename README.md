# BIDMC Respiratory-Rate Reference-Path Comparison

**A fixed-window technical audit of stored monitor respiratory rate versus two-annotator-derived respiratory rate in the BIDMC PPG and Respiration Dataset v1.0.0.**

Version: **v1.0.0** (research frozen 2026-10-13). Author: Sho Miyazaki (Independent researcher, Tokyo, Japan). Status: completed technical record, **not** a clinical validation or novel respiratory-rate estimation algorithm.

## What is included

For all 53 BIDMC records, compare the 1 Hz stored monitor `RESP` mean with a reference respiratory rate obtained independently per annotator from annotated event intervals, then averaged only when both annotator rates are computable and differ by **at most 2 breaths/min** before rounding. Inputs originate from the same impedance-respiration signal; this comparison is **not** an independent physiological ground truth.

The prespecified primary grid is **32 s windows, 3 s stride**; sensitivity is **64 s, 6 s**. The monitor needs a complete finite sample window. Finite zero is retained (not imputed as missing). Event sample IDs are 1-based at 125 Hz; events are sorted and exact duplicates removed per annotator; windows are `[start, end)`. Reported MAE is **numerical discordance**, not prediction error or medical-device accuracy. Overlapping windows are not independent observations.

| Metric | Primary (32 s / 3 s) | Sensitivity (64 s / 6 s) |
|---|---:|---:|
| Windows | 7,950 | 3,710 |
| Eligible | 7,271 | 3,503 |
| Window-pooled MAE (breaths/min) | 0.926271466 | 0.715722128 |
| Signed bias, monitor minus reference | +0.261569063 | +0.326401189 |
| Equal-record MAE | 1.247876058 (53 records) | 0.743510452 (52 records) |

Primary eligible zero-all windows: **11**, all from record 13; MAE **21.705972527** breaths/min. All-zero monitor values do not establish apnoea, device fault, or physiological cause. The sensitivity cohort excludes record 13 through the strict monitor completeness criterion, so a lower sensitivity MAE cannot be interpreted as an improvement caused by window length.

## Files and source

- `N2_Session3_output/N2_Session3_Window_Results.csv`, `N2_Session5_output/N2_Session5_Window_Results.csv`: all windows, including excluded windows and reasons.
- `N2_Session4_output/N2_Record_Summary.csv`, `N2_Session5_output/N2_Session5_Record_Summary.csv`: record-level summaries.
- Corresponding JSON summaries, normalization ledger and synthetic controls in session folders.
- The six top-level `N2_Session*.py` scripts are preserved historical source copies (not rewritten).
- `N2_Source_and_Input_Manifest.json`: source, upstream checksum and input provenance record.

**Raw CSVs are not redistributed.** Download from https://physionet.org/content/bidmc/1.0.0/ or run `python N2_Session1_Fetch_Verify.py` in the release root. The fetch script downloads the 106 official CSVs and verifies their sizes and SHA-256 against the official checksum list; network access is required for this first step. Do not overwrite original results with unverified outputs.

## Reproduce

Requires Python 3.9+ (standard library only), network access for download, and sufficient space for derived outputs. From the **release root**:

```bash
python N2_Session1_Fetch_Verify.py
python N2_Session2_RR_Controls_FIXED.py --input-dir N2_inputs_v1.0.0 --output-dir N2_Session2_output
python N2_Session3_Primary_Windows.py --input-dir N2_inputs_v1.0.0 --session2-ledger N2_Session2_output/N2_Session2_Normalization_Ledger.json --output-dir N2_Session3_output
cp N2_Session3_output/N2_Session3_Window_Results.csv .
cp N2_Session3_output/N2_Session3_Controls.json .
python N2_Session4_Primary_Summary.py
python N2_Session5_Sensitivity.py --input-dir N2_inputs_v1.0.0 --output-dir N2_Session5_output
python N2_Session6_Reproducibility.py
```

**Important:** The archived Session 4 script reads `N2_Session3_Window_Results.csv` and `N2_Session3_Controls.json` **at repository root**, whereas Session 3 writes in `N2_Session3_output/`. Before running Session 4, copy both current Session 3 outputs to repository root (they are included in the release and initially match):

```bash
cp N2_Session3_output/N2_Session3_Window_Results.csv .
cp N2_Session3_output/N2_Session3_Controls.json .
```

On Windows PowerShell use `Copy-Item` for these two commands. The `cp` steps above are required only once; the repeated reminder documents why Session 4 uses root-level files. Session 6 creates `N2_Session6_output/replay_*`; use a fresh checkout or remove the prior replay output before re-running Session 6. Confirm hashes against `N2_Session3_output/N2_Session3_Controls.json`, `N2_Session5_output/N2_Session5_Sensitivity_Summary.json`, and `N2_Session6_Verification.json`.

**Known historical provenance limit:** The Session 5 summary has a null Session-4-hash field, because that optional link was not recorded. It is not evidence that the two analyses disagree. In the original study, Session 4 summary values were not used as sensitivity calculation inputs.

## Interpretation and contribution limits

Breath-interval RR, two-annotator screening and the 32/64 s windows are not newly invented here. The incremental contribution is the **traceable fixed-window table** for all 53 records: two reference paths, eligibility, missingness, finite zeros, exclusion reasons and reproducibility evidence. This does not establish a novel clinical discovery, a better reference method, a new RR model or any prediction-performance effect.

Record 23's source `Time` starts at -1 s and is not silently shifted. The historical pre-protocol audit's **40 raw-local-order-issue windows** and the later **47 normalization-candidate windows** use different definitions; neither means 40 or 47 valid rate estimates changed. The original paper's 97.5% acceptance figure and this audit's 92.94% refer to populations/rules not established as identical; this work **does not claim to exactly reproduce that figure**.

## Attribution, AI assistance, and citation

Upstream BIDMC v1.0.0: https://doi.org/10.13026/C2208R (ODC-By v1.0); Pimentel et al. (2017): https://doi.org/10.1109/TBME.2016.2613124. See `DATA_LICENSE_AND_ATTRIBUTION.md`. Author-created Python code is MIT licensed; this does not override data-origin rights.

ChatGPT assisted with protocol discussion, code generation, verification, independent audit and editorial organization. The named human author takes responsibility for final review, scientific interpretation, and distribution. AI is not an author.

**Citation:** Cite the version-specific Zenodo DOI once assigned; also cite the original BIDMC dataset and associated paper. This draft intentionally contains no unassigned or invented DOI.
