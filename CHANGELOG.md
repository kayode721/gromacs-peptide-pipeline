# Initial shareable package

Adapted from the author's existing `peptide_gmx_pipeline.sh` and input directory.

- Preserved the preparation → minimization → NVT → NPT → MD sequence and shared-preparation replicate workflow.
- Disabled postprocessing and plotting by default; DSSP remains optional in the same script.
- Added configurable terminal prompting and recorded standard input for `pdb2gmx`.
- Added output-directory protection and early validation for common configuration mistakes.
- Copied sibling topology includes when starting from a prepared structure.
- Used a matching protein-only structure and an explicit selection for optional DSSP.
- Added example configuration, documentation, orchestration tests, and GitHub Actions checks.

The bundled MDP examples retain the author's original parameters. Packaging is not a scientific reparameterization of those inputs.
