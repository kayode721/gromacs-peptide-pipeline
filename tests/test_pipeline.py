"""Orchestration tests; no GROMACS calculation is performed here."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
FAKE_GMX = r'''#!/usr/bin/env python3
import json, os, sys
from pathlib import Path
args = sys.argv[1:]
command = args[0]
with open(os.environ['GMX_CALL_LOG'], 'a') as log:
    log.write(json.dumps({'args': args, 'cwd': os.getcwd()}) + '\n')
def option(key):
    return args[args.index(key) + 1] if key in args else None
def create(name, content='mock\n'):
    if name:
        Path(name).write_text(content)
if command == 'pdb2gmx':
    create(option('-o')); create(option('-p')); create(option('-i'))
    if os.environ.get('EXPECT_MENU_INPUT'):
        assert sys.stdin.read() == '1\n2\n'
elif command == 'grompp':
    create(option('-o'))
elif command == 'mdrun':
    prefix = option('-deffnm')
    for ext in ['gro', 'cpt', 'xtc', 'edr', 'log']:
        create(prefix + '.' + ext)
elif command in ('editconf', 'solvate', 'genion', 'trjconv'):
    create(option('-o'))
elif command == 'dssp' and '-h' not in args:
    create(option('-o'), 'HHH~~~\nHH~~~~\n')
'''

class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='gmx pipeline test ')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.gmx = self.base / 'fake-gmx'
        self.gmx.write_text(FAKE_GMX)
        self.gmx.chmod(0o755)
        self.pdb = self.base / 'input.pdb'
        self.pdb.write_text('END\n')
        self.log = self.base / 'calls.jsonl'
        self.env = os.environ.copy()
        # Keep user pipeline variables from affecting the isolated test run.
        for key in ['PREFIX', 'WORK_DIR', 'WORK_PARENT', 'MDP_DIR', 'START_AT',
                    'REPLICAS', 'NVT_GEN_SEED', 'CHECK_ONLY', 'SEEDS', 'GMXLIB']:
            self.env.pop(key, None)
        self.env['GMX_CALL_LOG'] = str(self.log)

    def run_pipeline(self, extra='', check=False):
        config = self.base / 'config.yaml'
        config.write_text(f'''prefix: demo
input_pdb: "{self.pdb}"
mdp_dir: "{ROOT / 'examples/mdp'}"
work_parent: "{self.base / 'runs'}"
gmx: "{self.gmx}"
''' + extra)
        return subprocess.run(['bash', str(ROOT / 'peptide_gmx_pipeline.sh'),
                               '--config', str(config)] + (['--check'] if check else []),
                              cwd=self.base, env=self.env, text=True, capture_output=True)

    def calls(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()] if self.log.exists() else []

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_default_ends_with_production(self):
        self.assert_success(self.run_pipeline())
        commands = [c['args'][0] for c in self.calls()]
        self.assertEqual(commands, ['pdb2gmx', 'editconf', 'solvate', 'grompp', 'genion',
                                    'grompp', 'mdrun', 'grompp', 'mdrun', 'grompp',
                                    'mdrun', 'grompp', 'mdrun'])
        self.assertTrue((self.base / 'runs/demo/demo.xtc').exists())
        self.assertFalse((self.base / 'runs/demo/demo_noPBC.xtc').exists())
        self.assertIn('-ter', self.calls()[0]['args'])

    def test_optional_dssp_uses_matching_structure(self):
        self.assert_success(self.run_pipeline('postprocess: yes\ndssp: yes\n'))
        call = [c['args'] for c in self.calls() if c['args'][0] == 'dssp' and '-h' not in c['args']][0]
        self.assertEqual(call[call.index('-s') + 1], 'demo_noPBC.gro')
        self.assertIn('-sel', call)
        self.assertTrue((self.base / 'runs/demo/demo.dat').exists())

    def test_replicates_share_preparation_and_have_distinct_seeds(self):
        self.assert_success(self.run_pipeline('replicas: 2\n'))
        self.assertEqual(sum(c['args'][0] == 'pdb2gmx' for c in self.calls()), 1)
        for number, seed in [(1, 107919), (2, 115838)]:
            folder = self.base / f'runs/demo_rep{number}'
            self.assertTrue((folder / f'demo_rep{number}.xtc').exists())
            self.assertRegex((folder / 'input/nvt.mdp').read_text(), rf'gen_seed\s*=\s*{seed}')

    def test_recorded_termini_propagate_into_shared_preparation(self):
        answers = self.base / 'answers.txt'
        answers.write_text('1\n2\n')
        self.env['EXPECT_MENU_INPUT'] = '1'
        self.assert_success(self.run_pipeline(f'replicas: 2\npdb2gmx_input: "{answers}"\n'))

    def test_termini_no_omits_ter(self):
        self.assert_success(self.run_pipeline('termini: no\n'))
        self.assertNotIn('-ter', self.calls()[0]['args'])

    def test_check_does_not_launch_gromacs(self):
        self.assert_success(self.run_pipeline(check=True))
        self.assertEqual(self.calls(), [])
        self.assertFalse((self.base / 'runs').exists())

    def test_invalid_dssp_combination_fails_early(self):
        result = self.run_pipeline('dssp: yes\npostprocess: no\n')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('dssp requires', result.stderr)
        self.assertEqual(self.calls(), [])

    def test_existing_outputs_are_not_overwritten(self):
        self.assert_success(self.run_pipeline())
        count = len(self.calls())
        result = self.run_pipeline()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('not an empty directory', result.stderr)
        self.assertEqual(len(self.calls()), count)

    def test_stop_after_md_skips_analysis(self):
        self.assert_success(self.run_pipeline('stop_after: md\npostprocess: yes\n'))
        self.assertNotIn('trjconv', [c['args'][0] for c in self.calls()])

    def test_prepared_topology_copies_sibling_includes(self):
        top = self.base / 'topol.top'
        gro = self.base / 'prepared.gro'
        top.write_text('#include "chain.itp"\n')
        gro.write_text('mock\n')
        (self.base / 'chain.itp').write_text('mock\n')
        self.assert_success(self.run_pipeline(f'start_at: editconf\ntopology: "{top}"\nprocessed_gro: "{gro}"\n'))
        self.assertNotIn('pdb2gmx', [c['args'][0] for c in self.calls()])
        self.assertTrue((self.base / 'runs/demo/chain.itp').exists())

if __name__ == '__main__':
    unittest.main()
