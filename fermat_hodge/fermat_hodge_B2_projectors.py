
# fermat_hodge_B2_projectors.sage
# (c) 2025. Vanishing-cycles + character-projector route (Algorithm B2)
# SageMath >= 9.0
#
# Idea: Use the standard Pham/Brieskorn picture. The vanishing homology of
#   f = x_0^m + ... + x_{n+1}^m
# in degree n is isomorphic (for projective smooth X_m^n) to the primitive middle
# homology/cohomology. A Z-basis can be indexed by tuples u in {1,...,m-1}^{n+2};
# the intersection form is the tensor product of the A_{m-1} form (Cartan matrix).
# The diagonal symmetry G acts by characters χ_u(g) = ζ_m^{⟨u, t(g)⟩}.
# Projecting to the Hodge (p,p) part amounts to selecting exactly the tuples u
# in Shioda's set B_m^n.
#
# This file provides:
#   - list_B_tuples(m,n): enumerate the tuples u ∈ {1..m-1}^{n+2} in B_m^n
#   - pham_A1d(m): the 1D A_{m-1} Cartan form
#   - primitive_gram_vanishing(m,n): Gram matrix on the B_m^n-indexed primitive lattice
#   - smith_invariants_from_matrix(A): SNF invariants
#   - full_invariants_from_primitive(invariants, m): append m for h^p
#   - compare_with_linear_route_fourfold(m): sanity check for n=4 (requires fermat_hodge_lattice.sage)
#
from sage.all import *
from itertools import product

def _units_mod(m):
    return [t for t in range(1, m) if gcd(t, m) == 1]

def list_B_tuples(m, n):
    """
    Enumerate Shioda's B_m^n as tuples u = (u_0,...,u_{n+1}) with entries in {1,...,m-1}.
    Conditions:
      (i) sum u_i ≡ 0 (mod m)
      (ii) sum <t u_i> = (p+1) m  for all t ∈ (Z/mZ)^*, where <a> ∈ {1,...,m-1}.
    Returns a list of tuples.
    """
    if n % 2 != 0:
        raise ValueError("n must be even (n = 2p).")
    p = n // 2
    k = n + 2
    units = _units_mod(m)
    target = (p + 1) * m

    B = []
    for u in product(range(1, m), repeat=k):
        if sum(u) % m != 0:
            continue
        good = True
        for t in units:
            s = 0
            for ui in u:
                w = (t*ui) % m
                if w == 0:
                    # (should not happen since ui in 1..m-1 and gcd(t,m)=1)
                    w = m
                s += w
            if s != target:
                good = False
                break
        if good:
            B.append(u)
    return B

def pham_A1d(m):
    """
    The 1D vanishing-cycle intersection form for A_{m-1} as an (m-1)x(m-1) matrix.
    """
    M = Matrix(ZZ, m-1, m-1)
    for i in range(m-1):
        M[i,i] = 2
        if i+1 < m-1:
            M[i,i+1] = -1
            M[i+1,i] = -1
    return M

def _tensor_entry(u, v, A1d):
    """
    Entry of the tensor-product Gram on basis vectors indexed by u,v.
    u,v: tuples in {1,...,m-1}^{k}
    A1d: (m-1)x(m-1) Cartan matrix for A_{m-1}
    """
    prod = ZZ(1)
    for a, b in zip(u, v):
        prod *= A1d[a-1, b-1]
        if prod == 0:
            return ZZ(0)
    return ZZ(prod)

def primitive_gram_vanishing(m, n, return_B=False):
    """
    Build the primitive Gram matrix on the Hodge (p,p) sublattice using the vanishing-cycle picture.
    Steps:
      1) Enumerate B = B_m^n (tuples u).
      2) Build the Gram matrix G with entries G_ij = Π_{ℓ} A1d[u_i[ℓ], u_j[ℓ]],
         where A1d is the A_{m-1} form.
    Returns the integer matrix G. If return_B=True, also returns the list B.
    """
    if n % 2 != 0:
        raise ValueError("n must be even (n = 2p).")
    B = list_B_tuples(m, n)
    A1d = pham_A1d(m)
    r = len(B)
    G = Matrix(ZZ, r, r)
    for i in range(r):
        ui = B[i]
        G[i,i] = _tensor_entry(ui, ui, A1d)
        for j in range(i+1, r):
            val = _tensor_entry(ui, B[j], A1d)
            G[i,j] = G[j,i] = val
    return (G, B) if return_B else G

def smith_invariants_from_matrix(A):
    """
    Smith normal form invariant factors (nonzero diagonal entries).
    """
    U, S, V = A.smith_form()
    diag = [ZZ(S[i,i]) for i in range(min(S.nrows(), S.ncols()))]
    return [d for d in diag if d != 0]

def full_invariants_from_primitive(primitive_invariants, m):
    """
    From primitive invariants to full Hodge lattice: append <h^p,h^p>=m as a 1x1 block.
    (Orthogonal direct sum since primitive is orthogonal to the hyperplane power.)
    """
    inv = list(primitive_invariants) + [ZZ(m)]
    inv.sort()
    return inv

############################################################
# Optional comparison with the linear route (n=4)          #
############################################################

def compare_with_linear_route_fourfold(m, max_matchings=None, target_rank=None, seed=0, verbose=True):
    """
    For n=4:
      - Compute primitive Gram from vanishing cycles (B2), get SNF invariants.
      - Compute primitive Gram from linear cycles (B1) on a greedy full-rank subset, get SNF.
      - Compare the invariant lists (after appending m to both for the full lattice).
    Requires 'fermat_hodge_lattice.sage' to be loaded.
    """
    # B2 side
    G = primitive_gram_vanishing(m, 4)
    prim_inv_B2 = smith_invariants_from_matrix(G)
    full_inv_B2 = full_invariants_from_primitive(prim_inv_B2, m)

    # B1 side (uses the other module)
    if 'linear_cycles_SNF_fourfold' not in globals():
        raise RuntimeError("Please load('fermat_hodge_lattice.sage') before calling this comparison.")
    res = linear_cycles_SNF_fourfold(m, max_matchings=max_matchings, target_rank=target_rank, seed=seed)
    prim_inv_B1 = res['primitive_invariants']
    full_inv_B1  = res['full_invariants']

    ok = (prim_inv_B1 == prim_inv_B2) and (full_inv_B1 == full_inv_B2)
    if verbose:
        print(f"[compare m={m}] primitive invariants equal? {prim_inv_B1 == prim_inv_B2}")
        print(f"[compare m={m}] full invariants equal?      {full_inv_B1  == full_inv_B2}")
    return dict(ok=ok, prim_B1=prim_inv_B1, prim_B2=prim_inv_B2, full_B1=full_inv_B1, full_B2=full_inv_B2)
