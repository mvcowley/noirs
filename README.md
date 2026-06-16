# Noirs

## Noisy Immune Receptor Repertoire Simulator

Noirs is a Rust tool for generating synthetic adaptive immune receptor
sequencing (AIRR) data with realistic amplification and sequencing noise. The
simulator models the full process from repertoire sampling through PCR
amplification and sequencing, producing datasets that can be used to benchmark
analysis pipelines, error-correction methods, and repertoire reconstruction
algorithms.

Noirs main novelty stems from a "lazy" approach to the PCR simulation step where
a sparse selection of the tree branches are simulated.

> **Note**
>
> This repo is very much a work in progress!

## Features

- Simulate immune receptor abundance distributions using a Zipfian sampling
  model
- Model PCR amplification across multiple cycles
- Introduce PCR-induced mutations during amplification
- Simulate sequencing errors
- Generate per-UMI observation trees
- Export results as NumPy (`.npy`) arrays for downstream analysis

## Simulation Workflow

1. **Sampling**
   - Receptor abundances are drawn from a Zipf distribution.
   - Each simulated tree originates from a unique molecular identifier (UMI).

2. **PCR Amplification**
   - Molecules are amplified through a configurable number of PCR cycles.
   - Per-cycle amplification efficiencies can be specified.

3. **Sequencing**
   - Amplified molecules are sequenced with a configurable sequencing error
     rate.
   - Sequencing errors are modeled independently from PCR errors.

4. **Output**
   - Each simulated UMI is written as a NumPy array containing:
     - PCR lineage information
     - PCR mutation counts
     - Sequencing error counts

## Installation

### Prerequisites

- Rust (edition 2021)
- Cargo

### Build

```bash
git clone https://github.com/mvcowley/noirs.git
cd noirs
cargo build --release
```

## Usage

Run the simulator:

```bash
cargo run --release
```

Simulation parameters are currently configured in `src/main.rs`.

Example configuration:

```rust
let total_observations = 1_000_000;
let max_observations = 1e3;
let exponent = 2.0;

let reaction = pcr::Reaction {
    sites: 12,
    efficiencies: vec![0.95; 30],
    errors: vec![0.0001; 30],
};

let sequencer = sequence::Sequencer {
    error: 0.005,
};
```

Generated output files are written to the configured output directory:

```text
out/
├── 1.npy
├── 2.npy
├── 3.npy
└── ...
```

## Output Format

Each output file contains observations associated with a single simulated UMI.

Columns include:

| Column Group      | Description                               |
| ----------------- | ----------------------------------------- |
| Observation Tree  | PCR lineage information                   |
| PCR Mutations     | Mutations introduced during amplification |
| Sequencing Errors | Errors introduced during sequencing       |

Arrays are stored in NumPy `.npy` format for easy integration with Python-based
analysis workflows.

## Project Structure

```text
src/
├── amplicon.rs
├── fastx.rs
├── noise.rs
├── parse.rs
├── pcr.rs
├── sample.rs
├── sequence.rs
├── zipf.rs
└── main.rs
```

Key modules:

- `sample` – receptor abundance sampling
- `pcr` – PCR amplification and mutation simulation
- `sequence` – sequencing error simulation
- `main` – simulation configuration and execution

## Reproducibility

The simulator uses a seeded ChaCha8 random number generator, allowing
experiments to be reproduced exactly by reusing the same seed.
