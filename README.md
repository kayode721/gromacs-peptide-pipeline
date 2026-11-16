# GROMACS Peptide Pipeline

A configurable Bash workflow that connects the GROMACS preparation and simulation stages, from a prepared peptide PDB to a production XTC trajectory. Supply your structure, force field, MDP files, and a small configuration file. The script passes outputs between stages and keeps each run in its own folder.

```text
PDB → pdb2gmx → box → solvation → ions → minimization → NVT → NPT → production MD
                                                                          │
                                                                     XTC + TPR
                                                                          │ optional
                                                      centered protein XTC + GRO → DSSP → plot
```

**One pipeline, optional analysis.** Production MD is the default endpoint. Enable trajectory processing or DSSP in the same configuration when needed.

## Requirements

- Bash and standard Unix tools (`awk`, `find`, `cp`, `tr`).
- GROMACS installed and available as `gmx`, or set `gmx:` to its executable path.
- A prepared peptide/protein structure and a compatible installed force field.
- Five MDP files: `ions.mdp`, `minim.mdp`, `nvt.mdp`, `npt.mdp`, and `md.mdp`.
- Optional: PyMOL for fetching/selecting a peptide from a PDB ID.
- Optional: a GROMACS build with `gmx dssp` for DSSP analysis; Python 3 and Matplotlib for PNG plots.

The orchestration is based on an existing SAA peptide workflow. This is a peptide/protein-in-water workflow; ligand parameterization, membranes, and arbitrary mixed systems require adaptation. A successful script execution does not establish equilibration or scientific validity.

## Quick start

Clone this repository, then work from its root directory:

```bash
git clone https://github.com/kayode721/gromacs-peptide-pipeline.git
cd gromacs-peptide-pipeline
```

```bash
mkdir -p inputs
cp /path/to/your/prepared-peptide.pdb inputs/peptide.pdb
cp peptide_pipeline.example.yaml config.yaml
# Edit config.yaml: structure, force field, MDP directory, run length, termini, etc.
bash peptide_gmx_pipeline.sh --config config.yaml --check
bash peptide_gmx_pipeline.sh --config config.yaml
```

`--check` validates required paths and prints a plan. It does not run GROMACS, validate topology chemistry, or establish force-field compatibility. Actual runs stop on command failure and refuse to reuse a nonempty output directory.

The example configuration uses `amber99sb-ildn`, the original script default. The author's current SAA configuration uses `charmm36-feb2026_cgenff-5.0`. To use an external force field, set `force_field` to its name (without `.ff`) and `gmxlib` to the directory containing that `.ff` folder. Force-field files and research structures are not bundled.

### MDP examples

`examples/mdp/` contains the author's input templates: 100 ps NVT, 2 ns NPT, and 300 ns production with a 2 fs timestep. Production saves protein coordinates every 10 ps (`nstxout-compressed = 5000`). These files are reference inputs, not universal recommendations. Review nonbonded cutoffs/modifiers, dispersion corrections, restraints, coupling, timestep, and run lengths against your selected force field and system before use. In particular, switching from AMBER to CHARMM does not make an existing MDP suitable automatically.

Use `mdp_dir` to point to your own validated set. The pipeline copies that directory into each run. Equilibration templates still write coordinates and velocities to TRR files; production XTC size is only part of the total storage requirement.

## Terminal chemistry and unattended runs

By default, `termini: yes` passes `-ter` to `pdb2gmx`, which asks you to select termini. For the author's CHARMM36 amidated SAA peptide, use a prepared PDB ending at ALA27 (no separate `NH2`/`CT2` residue), then select the appropriate N-terminal patch and the C-terminal `CT2` amide patch. Patch names and menu numbers depend on the force field.

For a batch job, first inspect the menus for the exact input, force field, and GROMACS version, then save the required numeric responses, one per line, in a text file:

```yaml
termini: yes
pdb2gmx_input: ./inputs/termini.txt
```

The file supplies standard input to `pdb2gmx`. Account for every prompt, including additional chains. Do not reuse menu numbers blindly. `termini: no` uses force-field defaults; it does not mean a neutral terminus or the `None` patch. Other structure-dependent prompts may still occur.

## Outputs

With `prefix: peptide` and `work_parent: ./runs`, files are written under `runs/peptide/`:

- `peptide.xtc`: production coordinates for the group specified in `md.mdp`.
- `peptide.tpr`: production run input/topology.
- `peptide.gro`, `.cpt`, `.edr`, `.log`: final structure, restart checkpoint, energies, and run log.
- Preparation and equilibration files, topology includes, and a copied `input/` directory.

Keep the TPR and topology with the trajectory for downstream analysis; retain checkpoints for restarting. Protein-only XTC files contain fewer atoms than a full-system TPR, so downstream tools must use a matching selection/structure.

## Optional trajectory processing and DSSP

Edit the same config to enable centering and protein extraction:

```yaml
postprocess: yes
dssp: no
```

This creates `*_noPBC.xtc` and a matching `*_noPBC.gro` using `trjconv -pbc mol -center`. It makes molecules whole and centers the selected protein; it does not guarantee an unwrapped diffusion trajectory or reconstruct every multimeric assembly.

To include DSSP in that same run:

```yaml
postprocess: yes
dssp: yes
dssp_plot: no
```

DSSP uses the matching protein-only GRO/XTC pair and produces `PREFIX.dat`. This workflow uses the built-in `gmx dssp` command, not the legacy `gmx do_dssp` interface.

For optional plots:

```bash
python3 -m pip install -r requirements-analysis.txt
```

For headless compute nodes, set `export MPLBACKEND=Agg`. If your home cache is not writable, point `MPLCONFIGDIR` and `XDG_CACHE_HOME` to writable directories.

```yaml
dssp_plot: yes
python: python3
dssp_duration_ns: 300
dssp_subsample: 10
```

Set `dssp_duration_ns` to the actual duration represented by the trajectory. The plotting script distributes the supplied duration over the saved frames; it does not read XTC timestamps. Do not use this assumption for irregularly sampled or concatenated trajectories without checking the time axis. The supplied plot helper must be in the copied MDP directory at `input/plot_dssp_states.py`, or set `dssp_plot_script` to its path.

These options run analysis after MD. They are not an analysis-only resume mode; for an existing trajectory, run the downstream commands separately rather than rerunning production.

## Replicate simulations

```yaml
replicas: 3
start_index: 1
seeds: 107919,115838,123757
```

Preparation and minimization run once in the base folder. Each replicate then starts from the shared minimized system, receives copied inputs and its own NVT velocity seed, and runs NVT → NPT → MD sequentially. Seeds above are examples; provide one distinct positive seed per replicate. Without an explicit list, the script generates a deterministic sequence starting at 107919. Replicates are sequential, not parallel scheduler jobs.

## Start from a prepared topology

```yaml
start_at: editconf
processed_gro: /path/to/peptide_processed.gro
topology: /path/to/topol.top
# posre_itp: /path/to/posre.itp
```

The script copies the topology and sibling `.itp` files, then starts box setup. Keep external force fields accessible via `gmxlib`. Nested topology include directories are not copied automatically; arrange their include paths before running.

## Other configuration details

- Paths resolve relative to the launch directory, not the YAML file.
- This is a flat `key: value` parser, not a general YAML parser. Avoid nested structures, multiline values, and `#` inside values.
- Environment variables can configure a run; keys in a supplied config take precedence.
- `gpu_mode`: `auto`, `cpu`, `gpu`, or `gpu-pme`. Hardware support comes from your GROMACS build.
- `ion_conc: none` neutralizes only; a number such as `0.15` requests that salt concentration in mol/L in addition to neutralization.
- `stop_after`: `prep`, `em`, `nvt`, `npt`, `md`, or `none`. `md` stops before optional analysis; `none` runs every enabled stage.
- `maxwarn: 0` is the simulation preprocessing default. `ion_maxwarn: 1` retains the original workflow's allowance for a net-charge warning before neutralization. Inspect warnings; a count does not distinguish their causes.
- Stage stopping is not automatic checkpoint resume. Resume interrupted production using GROMACS checkpoint facilities after checking the existing run.
- Inspect minimization convergence, density/temperature/pressure behavior, and GROMACS warnings before relying on the production results. The pipeline does not decide whether equilibration is sufficient.

## Validation

```bash
bash -n peptide_gmx_pipeline.sh
python3 -m unittest discover -s tests -v
```

The automated tests use a fake GROMACS executable to check command ordering, default output behavior, opt-in analysis, input failures, and replicate handling. They test orchestration, not molecular physics. See `VALIDATION.md` for checks performed on this packaged version.

## Attribution and reuse

This repository packages the author's Bash pipeline and DSSP plot helper. Cite GROMACS, the selected force field, and the analysis methods as appropriate in scientific work. The pipeline and plot helper are distributed under the [MIT License](LICENSE). External force fields, GROMACS, and other dependencies retain their own licenses.
