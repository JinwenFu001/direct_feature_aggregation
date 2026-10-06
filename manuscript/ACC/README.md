# JASA Author Contributions Checklist

These materials accompany **A direct approach to tree-guided feature aggregation for high-dimensional regression** by Jinwen Fu, Aaron J. Molstad and Hui Zou.

- `acc_form_draft.Rmd`: editable checklist source, updated 5 October 2026.
- `JASA_ACC_draft.docx`: editable Word version of the same checklist.
- `JASA_ACC_draft.pdf`: PDF review copy.
- `historical_environment.txt`: environment information extracted from a previously saved analysis result, not a new reproduction run.

The planned upload contains a top-level `README.md` and four directories: `code/`, `data/`, `manuscript/` and `output/`. This ACC belongs in `manuscript/ACC/`. The top-level README documents the complete directory structure, software requirements, script-to-figure/table mapping, execution order and use of saved results.

The current package contains 200 saved results for each of experiments 1-5, 400 for experiment 6, and 200 for each of the eleven real-data analyses. Figure 3-5 plotting scripts and the prepared figures are present. The Figure 3 plotting script currently writes `fig2joint.pdf`, while the prepared figure is named `Figure 3.pdf`; the checklist records this remaining filename inconsistency.

Real-data inputs are directly under `data/`, with eleven R runners and shell launchers under `code/run code/real data study/`. `Table 1.R` summarizes their saved results; its continuous-outcome aggregation follows the original notebook's row expansion and its signal-to-noise column is not computed. `sinha_2016_binary/Figure 6.R` performs a separate refit using the Sinha input and `sinha_2016_genera.tsv`. Data access/preprocessing documentation, dictionaries, environment validation, numerical agreement and measured runtime remain incomplete. Public availability and full reproduction are not certified by this update.

Earlier working directories, local tools and internal preparation audits are outside the upload. The ACC PDF contains the reader-facing description and limitations; internal author tasks and historical audit findings are kept in `jasa_acc_draft/jasa_minimum_requirements_review.txt` outside the release package. That internal checklist is not a submission attachment. Copies of the ACC in earlier preparation locations are kept synchronized for convenience; use `manuscript/ACC/` for subsequent editing and packaging. This update does not publish the repository.

Form guidance: [JASA ACC form and instructions](https://jasa-acs.github.io/repro-guide/pages/acc.html) and [JASA reproducibility guide](https://jasa-acs.github.io/repro-guide/).
