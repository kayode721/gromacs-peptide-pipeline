#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash peptide_gmx_pipeline.sh --config peptide_pipeline.example.yaml
  bash peptide_gmx_pipeline.sh --config peptide_pipeline.example.yaml --check

You can also configure values with environment variables (YAML values take precedence):
  PREFIX=saa1-27nh2 INPUT_PDB=/path/to/saa1-27nh2.pdb MDP_DIR=/path/to/input \
    bash peptide_gmx_pipeline.sh
  START_AT=editconf PROCESSED_GRO=/path/to/saa1-13_processed.gro TOPOLOGY=/path/to/topol.top \
    STOP_AFTER=npt bash peptide_gmx_pipeline.sh --config peptide_pipeline.example.yaml

Required:
  PREFIX       Output prefix, e.g. saa1-27nh2.
  MDP_DIR      Source input directory containing ions.mdp, minim.mdp, nvt.mdp, npt.mdp, md.mdp,
               and optionally plot_dssp_states.py.
  INPUT_PDB    Prepared peptide PDB to simulate.
    or
  PDB_ID + CHAIN + RESI
               Fetch a structure in PyMOL and save a peptide selection first.
    or
  START_AT=editconf + PROCESSED_GRO + TOPOLOGY
               Start from an existing pdb2gmx output and run box through production MD.

Common options:
  --config FILE    Flat YAML file with keys like prefix: saa1-27nh2.
  --check          Validate paths and print the run plan without running.
  --help          Show this help.

Force-field location (optional YAML key):
  gmxlib: .  # Parent directory containing the <force_field>.ff folder.
  Relative paths resolve from the launch directory, before changing folders.
  Use an absolute path to launch from anywhere. Sets GMXLIB for this run.

Stops:
  STOP_AFTER=prep,em,nvt,npt,md,none

Replicates:
  Set replicas: 3 to run PREFIX_rep1, PREFIX_rep2, and PREFIX_rep3.
  Each replicate gets preparation files, a copied input/ folder, and a different NVT gen_seed.

Terminal selection:
  termini: yes enables pdb2gmx -ter (default); choose termini interactively.
  pdb2gmx_input: /path/to/answers.txt supplies recorded menu selections for batch runs.
  termini: no uses force-field defaults; it does not select the None patches.

Postprocessing (all disabled by default):
  postprocess: yes creates *_noPBC.xtc and a matching *_noPBC.gro for viewing.
  dssp: yes runs gmx dssp on *_noPBC.xtc.
  dssp_plot: yes runs input/plot_dssp_states.py on the DSSP .dat output.

By default the script creates a new run folder named PREFIX, copies MDP_DIR into
that folder as input/, and writes all GROMACS outputs inside the run folder.
USAGE
}

trim() {
  local s="$1"
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "$s"
}

strip_quotes() {
  local s
  s="$(trim "$1")"
  if [[ "$s" == \"*\" && "$s" == *\" ]]; then
    s="${s:1:${#s}-2}"
  elif [[ "$s" == \'*\' && "$s" == *\' ]]; then
    s="${s:1:${#s}-2}"
  fi
  printf '%s' "$s"
}

load_simple_yaml() {
  local file="$1"
  local line key value

  if [[ ! -f "$file" ]]; then
    echo "ERROR: config not found: $file" >&2
    exit 1
  fi

  while IFS= read -r line || [[ -n "$line" ]]; do
    line="${line%%#*}"
    line="$(trim "$line")"
    [[ -z "$line" ]] && continue

    if [[ "$line" != *:* ]]; then
      echo "ERROR: config line is not key: value syntax: $line" >&2
      exit 1
    fi

    key="$(trim "${line%%:*}")"
    value="$(strip_quotes "${line#*:}")"

    case "$key" in
      prefix) PREFIX="$value" ;;
      start_at) START_AT="$value" ;;
      input_pdb) INPUT_PDB="$value" ;;
      processed_gro) PROCESSED_GRO="$value" ;;
      topology) TOPOLOGY="$value" ;;
      posre_itp) POSRE_ITP="$value" ;;
      pdb_id) PDB_ID="$value" ;;
      chain) CHAIN="$value" ;;
      resi) RESI="$value" ;;
      work_parent) WORK_PARENT="$value" ;;
      work_dir) WORK_DIR="$value" ;;
      mdp_dir) MDP_DIR="$value" ;;
      copy_input) COPY_INPUT="$value" ;;
      force_field) FORCE_FIELD="$value" ;;
      gmxlib) CONFIG_GMXLIB="$value" ;;
      water) WATER="$value" ;;
      termini) TERMINI="$value" ;;
      pdb2gmx_input) PDB2GMX_INPUT="$value" ;;
      box_type) BOX_TYPE="$value" ;;
      box_distance) BOX_DISTANCE="$value" ;;
      ion_conc) ION_CONC="$value" ;;
      gmx) GMX="$value" ;;
      pymol) PYMOL="$value" ;;
      gpu_mode) GPU_MODE="$value" ;;
      maxwarn) MAXWARN="$value" ;;
      ion_maxwarn) ION_MAXWARN="$value" ;;
      replicas) REPLICAS="$value" ;;
      start_index) START_INDEX="$value" ;;
      seeds) SEEDS="$value" ;;
      stop_after) STOP_AFTER="$value" ;;
      postprocess) POSTPROCESS="$value" ;;
      dssp) DSSP="$value" ;;
      dssp_plot) DSSP_PLOT="$value" ;;
      dssp_duration_ns) DSSP_DURATION_NS="$value" ;;
      dssp_subsample) DSSP_SUBSAMPLE="$value" ;;
      dssp_plot_script) DSSP_PLOT_SCRIPT="$value" ;;
      python) PYTHON="$value" ;;
      *)
        echo "WARNING: ignoring unknown config key: $key" >&2
        ;;
    esac
  done < "$file"
}

as_abs_path() {
  local path="$1"
  if [[ "$path" == /* ]]; then
    printf '%s' "$path"
  else
    printf '%s/%s' "$CALL_DIR" "$path"
  fi
}

is_yes() {
  local value
  value="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
  case "$value" in
    1|yes|true|on) return 0 ;;
    *) return 1 ;;
  esac
}

lower() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

mdp_set() {
  local file="$1"
  local key="$2"
  local value="$3"
  local tmp="${file}.tmp"

  awk -v key="$key" -v value="$value" '
    BEGIN { found = 0 }
    {
      line = $0
      stripped = line
      sub(/^[[:space:]]*/, "", stripped)
      if (stripped ~ "^" key "[[:space:]]*=") {
        comment = ""
        if (index(line, ";") > 0) {
          comment = substr(line, index(line, ";"))
        }
        printf "%-23s = %s", key, value
        if (comment != "") {
          printf "        %s", comment
        }
        printf "\n"
        found = 1
      } else {
        print line
      }
    }
    END {
      if (found == 0) {
        printf "%-23s = %s\n", key, value
      }
    }
  ' "$file" > "$tmp"

  mv "$tmp" "$file"
}

seed_for_replica() {
  local zero_based="$1"
  local seed

  if [[ -n "$SEEDS" ]]; then
    IFS=',' read -r -a seed_array <<< "$SEEDS"
    if [[ "$zero_based" -ge "${#seed_array[@]}" ]]; then
      echo "ERROR: seeds has fewer entries than replicas." >&2
      exit 1
    fi
    seed="$(trim "${seed_array[$zero_based]}")"
    if [[ -z "$seed" ]]; then
      echo "ERROR: empty seed entry in seeds list." >&2
      exit 1
    fi
    printf '%s' "$seed"
  else
    printf '%s' "$((100000 + (zero_based + 1) * 7919))"
  fi
}

maybe_stop() {
  local stage="$1"
  if [[ "$STOP_AFTER" == "$stage" ]]; then
    echo "Stopped after $stage as requested."
    exit 0
  fi
}

run_mdrun() {
  local deffnm="$1"
  case "$GPU_MODE" in
    auto)
      "$GMX" mdrun -deffnm "$deffnm" -v
      ;;
    cpu)
      "$GMX" mdrun -deffnm "$deffnm" -nb cpu -v
      ;;
    gpu)
      "$GMX" mdrun -deffnm "$deffnm" -nb gpu -v
      ;;
    gpu-pme)
      "$GMX" mdrun -deffnm "$deffnm" -nb gpu -pme gpu -v
      ;;
    *)
      echo "ERROR: gpu_mode must be auto, cpu, gpu, or gpu-pme." >&2
      exit 1
      ;;
  esac
}

CALL_DIR="$(pwd)"

PREFIX="${PREFIX:-}"
START_AT="${START_AT:-pdb2gmx}"
INPUT_PDB="${INPUT_PDB:-}"
PROCESSED_GRO="${PROCESSED_GRO:-}"
TOPOLOGY="${TOPOLOGY:-}"
POSRE_ITP="${POSRE_ITP:-}"
PDB_ID="${PDB_ID:-}"
CHAIN="${CHAIN:-A}"
RESI="${RESI:-}"
WORK_PARENT="${WORK_PARENT:-}"
WORK_DIR="${WORK_DIR:-}"
MDP_DIR="${MDP_DIR:-}"
COPY_INPUT="${COPY_INPUT:-yes}"
FORCE_FIELD="${FORCE_FIELD:-amber99sb-ildn}"
WATER="${WATER:-tip3p}"
TERMINI="${TERMINI:-yes}"
PDB2GMX_INPUT="${PDB2GMX_INPUT:-}"
BOX_TYPE="${BOX_TYPE:-cubic}"
BOX_DISTANCE="${BOX_DISTANCE:-1.0}"
ION_CONC="${ION_CONC:-0.15}"
GMX="${GMX:-gmx}"
PYMOL="${PYMOL:-pymol}"
GPU_MODE="${GPU_MODE:-auto}"
MAXWARN="${MAXWARN:-0}"
ION_MAXWARN="${ION_MAXWARN:-1}"
REPLICAS="${REPLICAS:-1}"
START_INDEX="${START_INDEX:-1}"
SEEDS="${SEEDS:-}"
NVT_GEN_SEED="${NVT_GEN_SEED:-}"
STOP_AFTER="${STOP_AFTER:-none}"
POSTPROCESS="${POSTPROCESS:-no}"
DSSP="${DSSP:-no}"
DSSP_PLOT="${DSSP_PLOT:-no}"
DSSP_DURATION_NS="${DSSP_DURATION_NS:-300}"
DSSP_SUBSAMPLE="${DSSP_SUBSAMPLE:-10}"
DSSP_PLOT_SCRIPT="${DSSP_PLOT_SCRIPT:-input/plot_dssp_states.py}"
PYTHON="${PYTHON:-python3}"
CONFIG_GMXLIB=""
CHECK_ONLY="${CHECK_ONLY:-no}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      [[ $# -ge 2 ]] || { echo "ERROR: --config requires a file." >&2; exit 1; }
      load_simple_yaml "$2"
      shift 2
      ;;
    --check)
      CHECK_ONLY="yes"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

# Resolve YAML's library path before any run-directory changes.
if [[ -n "$CONFIG_GMXLIB" ]]; then
  GMXLIB="$(as_abs_path "$CONFIG_GMXLIB")"
  if [[ ! -f "$GMXLIB/$FORCE_FIELD.ff/forcefield.itp" ]]; then
    echo "ERROR: force-field file not found: $GMXLIB/$FORCE_FIELD.ff/forcefield.itp" >&2
    echo "Set gmxlib to the parent directory containing $FORCE_FIELD.ff." >&2
    exit 1
  fi
  export GMXLIB
fi

if [[ -z "$PREFIX" ]]; then
  echo "ERROR: set prefix in the config or PREFIX in the environment." >&2
  exit 1
fi

if ! [[ "$PREFIX" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
  echo "ERROR: prefix must start with a letter or digit and contain only letters, digits, _, ., or -." >&2
  exit 1
fi
case "$STOP_AFTER" in prep|em|nvt|npt|md|none) ;; *) echo "ERROR: invalid stop_after: $STOP_AFTER" >&2; exit 1 ;; esac
case "$GPU_MODE" in auto|cpu|gpu|gpu-pme) ;; *) echo "ERROR: invalid gpu_mode: $GPU_MODE" >&2; exit 1 ;; esac
if [[ -n "$PDB2GMX_INPUT" ]]; then
  PDB2GMX_INPUT="$(as_abs_path "$PDB2GMX_INPUT")"
  [[ -f "$PDB2GMX_INPUT" ]] || { echo "ERROR: pdb2gmx_input not found: $PDB2GMX_INPUT" >&2; exit 1; }
fi

if [[ -z "$MDP_DIR" ]]; then
  echo "ERROR: set mdp_dir to the folder containing your .mdp files." >&2
  exit 1
fi

if ! [[ "$REPLICAS" =~ ^[0-9]+$ ]] || [[ "$REPLICAS" -lt 1 ]]; then
  echo "ERROR: replicas must be a positive integer." >&2
  exit 1
fi

if ! [[ "$START_INDEX" =~ ^[0-9]+$ ]] || [[ "$START_INDEX" -lt 1 ]]; then
  echo "ERROR: start_index must be a positive integer." >&2
  exit 1
fi

START_AT="$(lower "$START_AT")"
case "$START_AT" in
  pdb2gmx|pdb)
    START_AT="pdb2gmx"
    ;;
  editconf|processed|box)
    START_AT="editconf"
    ;;
  *)
    echo "ERROR: start_at must be pdb2gmx or editconf." >&2
    exit 1
    ;;
esac

if [[ "$START_AT" == "pdb2gmx" ]]; then
  if [[ -z "$INPUT_PDB" && -z "$PDB_ID" ]]; then
    echo "ERROR: set input_pdb, or set pdb_id with chain and resi." >&2
    exit 1
  fi

  if [[ -n "$PDB_ID" && -z "$RESI" ]]; then
    echo "ERROR: pdb_id mode requires resi, e.g. resi: 1-13." >&2
    exit 1
  fi
else
  if [[ -z "$PROCESSED_GRO" || -z "$TOPOLOGY" ]]; then
    echo "ERROR: start_at: editconf requires processed_gro and topology." >&2
    exit 1
  fi
fi

if [[ -z "$WORK_DIR" ]]; then
  if [[ -n "$WORK_PARENT" ]]; then
    WORK_PARENT="$(as_abs_path "$WORK_PARENT")"
    WORK_DIR="$WORK_PARENT/$PREFIX"
  else
    WORK_DIR="$CALL_DIR/$PREFIX"
  fi
fi

WORK_DIR="$(as_abs_path "$WORK_DIR")"
MDP_DIR="$(as_abs_path "$MDP_DIR")"

if [[ -n "$INPUT_PDB" ]]; then
  INPUT_PDB="$(as_abs_path "$INPUT_PDB")"
  if [[ ! -f "$INPUT_PDB" ]]; then
    echo "ERROR: input_pdb not found: $INPUT_PDB" >&2
    exit 1
  fi
fi

if [[ -n "$PROCESSED_GRO" ]]; then
  PROCESSED_GRO="$(as_abs_path "$PROCESSED_GRO")"
  if [[ ! -f "$PROCESSED_GRO" ]]; then
    echo "ERROR: processed_gro not found: $PROCESSED_GRO" >&2
    exit 1
  fi
fi

if [[ -n "$TOPOLOGY" ]]; then
  TOPOLOGY="$(as_abs_path "$TOPOLOGY")"
  if [[ ! -f "$TOPOLOGY" ]]; then
    echo "ERROR: topology not found: $TOPOLOGY" >&2
    exit 1
  fi
fi

if [[ -n "$POSRE_ITP" ]]; then
  POSRE_ITP="$(as_abs_path "$POSRE_ITP")"
  if [[ ! -f "$POSRE_ITP" ]]; then
    echo "ERROR: posre_itp not found: $POSRE_ITP" >&2
    exit 1
  fi
fi

for mdp in ions.mdp minim.mdp nvt.mdp npt.mdp md.mdp; do
  if [[ ! -f "$MDP_DIR/$mdp" ]]; then
    echo "ERROR: missing MDP file: $MDP_DIR/$mdp" >&2
    exit 1
  fi
done

if [[ "$REPLICAS" -gt 1 || -n "$NVT_GEN_SEED" ]] && ! is_yes "$COPY_INPUT"; then
  echo "ERROR: copy_input must be yes when running replicates with generated NVT seeds." >&2
  exit 1
fi

if is_yes "$DSSP" && ! is_yes "$POSTPROCESS"; then
  echo "ERROR: dssp requires postprocess: yes so the *_noPBC.xtc trajectory exists." >&2
  exit 1
fi

if is_yes "$CHECK_ONLY"; then
  echo "Run plan:"
  echo "  prefix: $PREFIX"
  echo "  start_at: $START_AT"
  if [[ "$START_AT" == "pdb2gmx" ]]; then
    echo "  input_pdb: ${INPUT_PDB:-<PyMOL fetch $PDB_ID chain $CHAIN resi $RESI>}"
  else
    echo "  processed_gro: $PROCESSED_GRO"
    echo "  topology: $TOPOLOGY"
    echo "  posre_itp: ${POSRE_ITP:-<copy from topology folder if present>}"
  fi
  echo "  work_dir: $WORK_DIR"
  echo "  source_input_dir: $MDP_DIR"
  if is_yes "$COPY_INPUT"; then
    echo "  run_input_dir: $WORK_DIR/input"
  else
    echo "  run_input_dir: <uses source input dir directly>"
  fi
  echo "  force_field: $FORCE_FIELD"
  echo "  GMXLIB: ${GMXLIB:-GROMACS default}"
  echo "  water: $WATER"
  echo "  termini: $TERMINI"
  echo "  pdb2gmx_input: ${PDB2GMX_INPUT:-<interactive/defaults>}"
  echo "  box: $BOX_TYPE, d=$BOX_DISTANCE nm"
  if [[ -n "$ION_CONC" && "$(lower "$ION_CONC")" != "none" ]]; then
    echo "  ions: $ION_CONC M, neutralized"
  else
    echo "  ions: neutralize only"
  fi
  echo "  gpu_mode: $GPU_MODE"
  echo "  replicas: $REPLICAS"
  echo "  start_index: $START_INDEX"
  echo "  stop_after: $STOP_AFTER"
  echo "  postprocess: $POSTPROCESS"
  echo "  dssp: $DSSP"
  if is_yes "$DSSP"; then
    echo "  dssp_plot: $DSSP_PLOT"
    echo "  python: $PYTHON"
    echo "  dssp_duration_ns: $DSSP_DURATION_NS"
    echo "  dssp_subsample: $DSSP_SUBSAMPLE"
    echo "  dssp_plot_script: $DSSP_PLOT_SCRIPT"
  fi
  if [[ "$REPLICAS" -gt 1 ]]; then
    if [[ -n "$WORK_PARENT" ]]; then
      replicate_parent="$WORK_PARENT"
    else
      replicate_parent="$(dirname "$WORK_DIR")"
    fi
    echo "  common prepared system:"
    echo "    $WORK_DIR"
    echo "    builds through neutralization and EM once; replicas start from ${PREFIX}_em.gro"
    echo "  replicate folders:"
    for ((i = 0; i < REPLICAS; i++)); do
      rep_num=$((START_INDEX + i))
      rep_prefix="${PREFIX}_rep${rep_num}"
      seed="$(seed_for_replica "$i")"
      echo "    ${rep_prefix}: seed $seed -> $replicate_parent/$rep_prefix"
    done
  elif [[ -n "$NVT_GEN_SEED" ]]; then
    echo "  nvt_gen_seed: $NVT_GEN_SEED"
  fi
  exit 0
fi

command -v "$GMX" >/dev/null 2>&1 || { echo "ERROR: GROMACS executable not found: $GMX" >&2; exit 1; }
if is_yes "$DSSP"; then
  "$GMX" dssp -h >/dev/null 2>&1 || { echo "ERROR: this GROMACS build does not provide gmx dssp." >&2; exit 1; }
fi
if is_yes "$DSSP" && is_yes "$DSSP_PLOT"; then
  command -v "$PYTHON" >/dev/null 2>&1 || { echo "ERROR: Python executable not found: $PYTHON" >&2; exit 1; }
fi
assert_empty_destination() {
  local destination="$1"
  if [[ -e "$destination" ]] && { [[ ! -d "$destination" ]] || [[ -n "$(ls -A "$destination")" ]]; }; then
    echo "ERROR: output destination is not an empty directory: $destination. Choose a new prefix/work_dir." >&2
    exit 1
  fi
}
assert_empty_destination "$WORK_DIR"
if [[ "$REPLICAS" -gt 1 ]]; then
  check_parent="${WORK_PARENT:-$(dirname "$WORK_DIR")}"
  for ((i = 0; i < REPLICAS; i++)); do
    assert_empty_destination "$check_parent/${PREFIX}_rep$((START_INDEX + i))"
    seed_for_replica "$i" >/dev/null
  done
fi

if [[ "$REPLICAS" -gt 1 ]]; then
  if [[ -n "$WORK_PARENT" ]]; then
    REPLICATE_PARENT="$WORK_PARENT"
  else
    REPLICATE_PARENT="$(dirname "$WORK_DIR")"
  fi

  BASE_PREFIX="$PREFIX"
  BASE_DIR="$WORK_DIR"
  case "$STOP_AFTER" in
    prep|em)
      COMMON_STOP_AFTER="$STOP_AFTER"
      ;;
    *)
      COMMON_STOP_AFTER="em"
      ;;
  esac

  mkdir -p "$REPLICATE_PARENT"

  echo
  echo "=== Building one shared processed/neutralized/minimized system ==="
  echo "Common system directory: $BASE_DIR"
  PREFIX="$BASE_PREFIX" \
  START_AT="$START_AT" \
  INPUT_PDB="$INPUT_PDB" \
  PROCESSED_GRO="$PROCESSED_GRO" \
  TOPOLOGY="$TOPOLOGY" \
  POSRE_ITP="$POSRE_ITP" \
  PDB_ID="$PDB_ID" \
  CHAIN="$CHAIN" \
  RESI="$RESI" \
  WORK_DIR="$BASE_DIR" \
  WORK_PARENT="" \
  MDP_DIR="$MDP_DIR" \
  COPY_INPUT="$COPY_INPUT" \
  FORCE_FIELD="$FORCE_FIELD" \
  WATER="$WATER" \
  TERMINI="$TERMINI" \
  PDB2GMX_INPUT="$PDB2GMX_INPUT" \
  BOX_TYPE="$BOX_TYPE" \
  BOX_DISTANCE="$BOX_DISTANCE" \
  ION_CONC="$ION_CONC" \
  GMX="$GMX" \
  PYMOL="$PYMOL" \
  GPU_MODE="$GPU_MODE" \
  MAXWARN="$MAXWARN" \
  ION_MAXWARN="$ION_MAXWARN" \
  REPLICAS=1 \
  NVT_GEN_SEED="" \
  STOP_AFTER="$COMMON_STOP_AFTER" \
  POSTPROCESS=no \
  DSSP=no \
  DSSP_PLOT=no \
  DSSP_DURATION_NS="$DSSP_DURATION_NS" \
  DSSP_SUBSAMPLE="$DSSP_SUBSAMPLE" \
  DSSP_PLOT_SCRIPT="$DSSP_PLOT_SCRIPT" \
  PYTHON="$PYTHON" \
  bash "${BASH_SOURCE[0]}"

  if [[ "$STOP_AFTER" == "prep" || "$STOP_AFTER" == "em" ]]; then
    echo
    echo "Stopped after shared $STOP_AFTER stage as requested."
    exit 0
  fi

  if [[ ! -f "$BASE_DIR/${BASE_PREFIX}_em.gro" ]]; then
    echo "ERROR: shared EM output not found: $BASE_DIR/${BASE_PREFIX}_em.gro" >&2
    exit 1
  fi

  if [[ ! -f "$BASE_DIR/topol.top" ]]; then
    echo "ERROR: shared topology not found: $BASE_DIR/topol.top" >&2
    exit 1
  fi

  for ((i = 0; i < REPLICAS; i++)); do
    rep_num=$((START_INDEX + i))
    rep_prefix="${PREFIX}_rep${rep_num}"
    rep_dir="$REPLICATE_PARENT/$rep_prefix"
    seed="$(seed_for_replica "$i")"

    echo
    echo "=== Launching $rep_prefix from shared ${BASE_PREFIX}_em.gro with NVT seed $seed ==="
    mkdir -p "$rep_dir/input"

    # Retain shared preparation outputs alongside each independent replicate.
    # Rename the prefix only; these files describe the same shared preparation.
    for prep_file in "$BASE_DIR/${BASE_PREFIX}"_{processed.gro,box.gro,solv.gro,ions.tpr,solv_ions.gro,em.*}; do
      [[ -f "$prep_file" ]] || continue
      prep_suffix="${prep_file##*/}"
      prep_suffix="${prep_suffix#"$BASE_PREFIX"}"
      cp "$prep_file" "$rep_dir/${rep_prefix}${prep_suffix}"
    done
    cp "$BASE_DIR/${BASE_PREFIX}_em.gro" "$rep_dir/${rep_prefix}_em.gro"
    cp "$BASE_DIR/topol.top" "$rep_dir/topol.top"
    find "$BASE_DIR" -maxdepth 1 -type f -name '*.itp' -exec cp {} "$rep_dir/" \;
    cp -R "$MDP_DIR"/. "$rep_dir/input/"

    run_mdp_dir="$rep_dir/input"
    mdp_set "$run_mdp_dir/nvt.mdp" "gen_vel" "yes"
    mdp_set "$run_mdp_dir/nvt.mdp" "gen_seed" "$seed"

    cd "$rep_dir"

    "$GMX" grompp \
      -f "$run_mdp_dir/nvt.mdp" \
      -c "${rep_prefix}_em.gro" \
      -r "${rep_prefix}_em.gro" \
      -p topol.top \
      -o "${rep_prefix}_nvt.tpr" \
      -maxwarn "$MAXWARN"
    run_mdrun "${rep_prefix}_nvt"
    if [[ "$STOP_AFTER" == "nvt" ]]; then
      echo "Stopped after nvt as requested for $rep_prefix."
      continue
    fi

    "$GMX" grompp \
      -f "$run_mdp_dir/npt.mdp" \
      -c "${rep_prefix}_nvt.gro" \
      -r "${rep_prefix}_nvt.gro" \
      -t "${rep_prefix}_nvt.cpt" \
      -p topol.top \
      -o "${rep_prefix}_npt.tpr" \
      -maxwarn "$MAXWARN"
    run_mdrun "${rep_prefix}_npt"
    if [[ "$STOP_AFTER" == "npt" ]]; then
      echo "Stopped after npt as requested for $rep_prefix."
      continue
    fi

    "$GMX" grompp \
      -f "$run_mdp_dir/md.mdp" \
      -c "${rep_prefix}_npt.gro" \
      -t "${rep_prefix}_npt.cpt" \
      -p topol.top \
      -o "${rep_prefix}.tpr" \
      -maxwarn "$MAXWARN"
    run_mdrun "$rep_prefix"
    if [[ "$STOP_AFTER" == "md" ]]; then
      echo "Stopped after md as requested for $rep_prefix."
      continue
    fi

    if is_yes "$POSTPROCESS"; then
      echo "Removing periodic boundary artifacts and centering $rep_prefix."
      printf 'Protein\nProtein\n' | "$GMX" trjconv \
        -s "${rep_prefix}.tpr" \
        -f "${rep_prefix}.xtc" \
        -o "${rep_prefix}_noPBC.xtc" \
        -pbc mol \
        -center
      # Use the identical GROMACS selection and frame processing for the structure.
      # This retains terminal caps included in Protein and preserves atom order.
      printf 'Protein\nProtein\n' | "$GMX" trjconv \
        -s "${rep_prefix}.tpr" \
        -f "${rep_prefix}.xtc" \
        -o "${rep_prefix}_noPBC.gro" \
        -pbc mol \
        -center \
        -dump 0
    fi

    if is_yes "$DSSP"; then
      echo "Running DSSP for $rep_prefix using ${rep_prefix}_noPBC.xtc."
      "$GMX" dssp \
        -s "${rep_prefix}_noPBC.gro" \
        -f "${rep_prefix}_noPBC.xtc" \
        -sel 'group "Protein"' \
        -o "${rep_prefix}.dat"

      if is_yes "$DSSP_PLOT"; then
        if [[ ! -f "$DSSP_PLOT_SCRIPT" ]]; then
          echo "ERROR: DSSP plot script not found in $rep_dir: $DSSP_PLOT_SCRIPT" >&2
          exit 1
        fi

        echo "Plotting DSSP secondary-structure map for $rep_prefix."
        "$PYTHON" "$DSSP_PLOT_SCRIPT" "${rep_prefix}.dat" \
          --duration-ns "$DSSP_DURATION_NS" \
          --subsample "$DSSP_SUBSAMPLE" \
          --title "${rep_prefix} DSSP secondary structure" \
          --out "${rep_prefix}_dssp_states_${DSSP_DURATION_NS}ns.png"
      fi
    fi

    echo "Done with $rep_prefix."
  done

  echo
  echo "All requested replicates finished."
  exit 0
fi

mkdir -p "$WORK_DIR"

RUN_MDP_DIR="$MDP_DIR"
if is_yes "$COPY_INPUT"; then
  mkdir -p "$WORK_DIR/input"
  cp -R "$MDP_DIR"/. "$WORK_DIR/input/"
  RUN_MDP_DIR="$WORK_DIR/input"
fi

if [[ -n "$NVT_GEN_SEED" ]]; then
  mdp_set "$RUN_MDP_DIR/nvt.mdp" "gen_vel" "yes"
  mdp_set "$RUN_MDP_DIR/nvt.mdp" "gen_seed" "$NVT_GEN_SEED"
fi

cd "$WORK_DIR"

echo "Prefix: $PREFIX"
echo "Work directory: $WORK_DIR"
echo "Input directory: $RUN_MDP_DIR"
echo "Force field: $FORCE_FIELD"

if [[ "$START_AT" == "pdb2gmx" ]]; then
  if [[ -n "$INPUT_PDB" ]]; then
    cp "$INPUT_PDB" "$PREFIX.pdb"
  else
    pml_file="prep_${PREFIX}.pml"
    {
      printf 'fetch %s, full, async=0\n' "$PDB_ID"
      printf 'remove solvent\n'
      printf 'create pep, full and chain %s and resi %s\n' "$CHAIN" "$RESI"
      printf 'save %s.pdb, pep\n' "$PREFIX"
      printf 'quit\n'
    } > "$pml_file"
    "$PYMOL" -cq "$pml_file"
  fi

  pdb2gmx_args=(
    -f "$PREFIX.pdb"
    -o "${PREFIX}_processed.gro"
    -p topol.top
    -i posre.itp
    -ff "$FORCE_FIELD"
    -water "$WATER"
    -ignh
  )
  if is_yes "$TERMINI"; then
    pdb2gmx_args+=(-ter)
  fi
  if [[ -n "$PDB2GMX_INPUT" ]]; then
    "$GMX" pdb2gmx "${pdb2gmx_args[@]}" < "$PDB2GMX_INPUT"
  else
    "$GMX" pdb2gmx "${pdb2gmx_args[@]}"
  fi
else
  cp "$PROCESSED_GRO" "${PREFIX}_processed.gro"
  cp "$TOPOLOGY" topol.top

  topology_dir="$(dirname "$TOPOLOGY")"
  find "$topology_dir" -maxdepth 1 -type f -name '*.itp' -exec cp {} "$WORK_DIR/" \;
  if [[ -n "$POSRE_ITP" ]]; then
    cp "$POSRE_ITP" posre.itp
  elif [[ -f "$topology_dir/posre.itp" ]]; then
    cp "$topology_dir/posre.itp" posre.itp
  fi
fi

"$GMX" editconf \
  -f "${PREFIX}_processed.gro" \
  -o "${PREFIX}_box.gro" \
  -c \
  -d "$BOX_DISTANCE" \
  -bt "$BOX_TYPE"

"$GMX" solvate \
  -cp "${PREFIX}_box.gro" \
  -cs spc216.gro \
  -o "${PREFIX}_solv.gro" \
  -p topol.top

"$GMX" grompp \
  -f "$RUN_MDP_DIR/ions.mdp" \
  -c "${PREFIX}_solv.gro" \
  -p topol.top \
  -o "${PREFIX}_ions.tpr" \
  -maxwarn "$ION_MAXWARN"

genion_args=(
  -s "${PREFIX}_ions.tpr"
  -o "${PREFIX}_solv_ions.gro"
  -p topol.top
  -pname NA
  -nname CL
  -neutral
)

if [[ -n "$ION_CONC" && "$(lower "$ION_CONC")" != "none" ]]; then
  genion_args+=(-conc "$ION_CONC")
fi

printf 'SOL\n' | "$GMX" genion "${genion_args[@]}"

maybe_stop prep

"$GMX" grompp \
  -f "$RUN_MDP_DIR/minim.mdp" \
  -c "${PREFIX}_solv_ions.gro" \
  -p topol.top \
  -o "${PREFIX}_em.tpr" \
  -maxwarn "$MAXWARN"
run_mdrun "${PREFIX}_em"
maybe_stop em

"$GMX" grompp \
  -f "$RUN_MDP_DIR/nvt.mdp" \
  -c "${PREFIX}_em.gro" \
  -r "${PREFIX}_em.gro" \
  -p topol.top \
  -o "${PREFIX}_nvt.tpr" \
  -maxwarn "$MAXWARN"
run_mdrun "${PREFIX}_nvt"
maybe_stop nvt

"$GMX" grompp \
  -f "$RUN_MDP_DIR/npt.mdp" \
  -c "${PREFIX}_nvt.gro" \
  -r "${PREFIX}_nvt.gro" \
  -t "${PREFIX}_nvt.cpt" \
  -p topol.top \
  -o "${PREFIX}_npt.tpr" \
  -maxwarn "$MAXWARN"
run_mdrun "${PREFIX}_npt"
maybe_stop npt

"$GMX" grompp \
  -f "$RUN_MDP_DIR/md.mdp" \
  -c "${PREFIX}_npt.gro" \
  -t "${PREFIX}_npt.cpt" \
  -p topol.top \
  -o "${PREFIX}.tpr" \
  -maxwarn "$MAXWARN"
run_mdrun "$PREFIX"
maybe_stop md

if is_yes "$POSTPROCESS"; then
  echo "Removing periodic boundary artifacts and centering $PREFIX."
  printf 'Protein\nProtein\n' | "$GMX" trjconv \
    -s "${PREFIX}.tpr" \
    -f "${PREFIX}.xtc" \
    -o "${PREFIX}_noPBC.xtc" \
    -pbc mol \
    -center
  # Use the identical GROMACS selection and frame processing for the structure.
  # This retains terminal caps included in Protein and preserves atom order.
  printf 'Protein\nProtein\n' | "$GMX" trjconv \
    -s "${PREFIX}.tpr" \
    -f "${PREFIX}.xtc" \
    -o "${PREFIX}_noPBC.gro" \
    -pbc mol \
    -center \
    -dump 0
fi

if is_yes "$DSSP"; then
  echo "Running DSSP for $PREFIX using ${PREFIX}_noPBC.xtc."
  "$GMX" dssp \
    -s "${PREFIX}_noPBC.gro" \
    -f "${PREFIX}_noPBC.xtc" \
    -sel 'group "Protein"' \
    -o "${PREFIX}.dat"

  if is_yes "$DSSP_PLOT"; then
    if [[ ! -f "$DSSP_PLOT_SCRIPT" ]]; then
      echo "ERROR: DSSP plot script not found in $WORK_DIR: $DSSP_PLOT_SCRIPT" >&2
      exit 1
    fi

    echo "Plotting DSSP secondary-structure map for $PREFIX."
    "$PYTHON" "$DSSP_PLOT_SCRIPT" "${PREFIX}.dat" \
      --duration-ns "$DSSP_DURATION_NS" \
      --subsample "$DSSP_SUBSAMPLE" \
      --title "${PREFIX} DSSP secondary structure" \
      --out "${PREFIX}_dssp_states_${DSSP_DURATION_NS}ns.png"
  fi
fi

echo "Done. Outputs are in $WORK_DIR."
