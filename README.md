# Fermat-Hodge

A SageMath toolkit for computing ranks and discriminants of Hodge lattices on Fermat hypersurfaces.

## Overview

This package provides computational tools for studying the Hodge lattices arising from Fermat hypersurfaces X_m^n, implementing algorithms based on:
- Shioda-Ran theory for computing ranks via character enumeration
- Linear cycles and intersection theory for discriminant computation
- Vanishing cycles and character projectors (Algorithm B2)

## Installation

### Prerequisites

- Python >= 3.7
- SageMath >= 9.0

### Install from source

```bash
git clone https://github.com/yourusername/fermat-hodge.git
cd fermat-hodge
pip install -e .
```

## Quick Start

### Using the SageMath modules directly

```sage
# In a SageMath session or notebook
load('fermat_hodge_lattice.sage')
load('fermat_hodge_B2_projectors.sage')

# Compute the rank of the Hodge lattice for X_5^4
m, n = 5, 4
rank = hodge_rank_via_B(m, n)
print(f"Rank of H^{{2,2}}_prim(X_{m}^{n}) = {rank}")

# Generate B-tuples
B = list_B_tuples(5, 2)
print(f"Number of B-tuples for X_5^2: {len(B)}")

# Compute Gram matrix via vanishing cycles
G = primitive_gram_vanishing(5, 2)
print(f"Gram matrix shape: {G.nrows()}x{G.ncols()}")
```

### Command Line Interface

The package includes a CLI for computing Hodge lattice invariants:

```bash
# Run with specific parameters
sage fermat_hodge_cli_runner.sage --m 5 --b1 --b2-m5

# Run comprehensive tests
sage test_suite.sage
```

## Package Structure

```
fermat-hodge/
├── fermat_hodge/
│   ├── __init__.py
│   ├── fermat_hodge_lattice.py       # Core rank and discriminant algorithms
│   ├── fermat_hodge_B2_projectors.py # Vanishing cycles approach
│   └── cli_runner.py                  # CLI interface
├── paper/                              # LaTeX documentation
│   ├── hodge_lattice_fermat_methods.tex
│   └── hodge_lattice_fermat_methods_with_appendix.tex
├── setup.py
└── README.md
```

## Main Functions

### Core Algorithms

- `hodge_rank_via_B(m, n)`: Compute rank using Shioda's character set B_m^n
- `hodge_discriminant_via_linear_cycles(m, n)`: Compute discriminant via linear cycles
- `hodge_invariants_linear_route(m, n)`: Get full Smith normal form invariants

### Vanishing Cycles Methods

- `list_B_tuples(m, n)`: Enumerate tuples in Shioda's set B_m^n
- `primitive_gram_vanishing(m, n)`: Gram matrix via vanishing cycles
- `smith_invariants_from_matrix(A)`: Compute Smith normal form invariants

## Mathematical Background

The package computes invariants of the primitive middle Hodge lattice H^{p,p}_prim(X_m^n) for Fermat hypersurfaces:

X_m^n: x_0^m + x_1^m + ... + x_{n+1}^m = 0

in P^{n+1}, where n = 2p is even.

## References

See the papers in the `paper/` directory for detailed mathematical exposition and proofs.

## License

[Specify your license here]

## Author

[Your Name]

## Contributing

Contributions are welcome! Please submit pull requests or open issues for bug reports and feature requests.