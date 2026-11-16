# Validation record

Checked on November 16, 2026.

## Automated checks

- Bash syntax check passed.
- 10 Python unittest orchestration tests passed using a fake GROMACS executable.
- Coverage includes default completion at production XTC, optional DSSP with a matching structure, shared preparation for replicates, distinct NVT seeds, terminal-answer propagation, disabled terminal prompting, check-only mode, invalid analysis settings, existing-output protection, and topology include copying.

## Real GROMACS execution check

A temporary SAA1–27 input was processed with GROMACS 2026.3 (Homebrew, CPU) and the author's external `charmm36-feb2026_cgenff-5.0` force field. The standalone residue-28 cap was removed from the temporary input and the CHARMM36 NH3+ and CT2 terminal patches were selected through a recorded answer file. The original research input was not changed or included in this repository.

The packaged script completed:

1. Topology generation with an amidated C-terminus.
2. Cubic box setup, solvation, and the ion preparation stage.
3. Energy minimization (maximum 500 steps).
4. NVT and NPT execution (20 steps each).
5. Production execution (20 steps at 2 fs), producing a protein-only XTC.
6. Protein centering/whole-molecule processing and a matching GRO snapshot.
7. Built-in GROMACS DSSP, analyzing three frames through 0.040 ps.

The solvated test system contained 25,026 atoms, including 432 protein atoms. This test confirms command and file compatibility for this example. It does not establish equilibration, convergence, scientifically appropriate production parameters, or successful completion of a 300 ns run. The GROMACS force field and test structure are external dependencies.

DSSP with a GRO structure emitted GROMACS' notice about inferred masses/radii. Other versions and molecular systems have not been certified. The optional plotting helper was also run separately on the three-frame DSSP result and produced a PNG and summary CSV. It used the Agg backend and temporary writable Matplotlib/font caches because the execution sandbox cannot write to the default user cache directories.
