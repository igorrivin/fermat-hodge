"""
Fermat-Hodge: Toolkit for ranks and discriminants of Hodge lattices on Fermat hypersurfaces.

This package provides computational tools for studying Hodge lattices on Fermat hypersurfaces
using SageMath. It implements algorithms for computing ranks and discriminants via different
approaches including Shioda-Ran theory and vanishing cycles.
"""

__version__ = "0.1.0"

# Note: This package contains SageMath code (.sage files converted to .py)
# To use these modules, you need to have SageMath installed and run them
# within a SageMath environment.

# The main modules are:
# - fermat_hodge_lattice: Core rank and discriminant algorithms  
# - fermat_hodge_B2_projectors: Vanishing cycles approach
# - fermat_hodge_cli_runner: Command-line interface
# - driver_fourfold_rank_logger: Helper for streaming computations

# To use the modules directly in SageMath:
# sage: load('fermat_hodge_lattice.sage')
# sage: load('fermat_hodge_B2_projectors.sage')