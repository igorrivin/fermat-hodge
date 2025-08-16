
# fermat_hodge_lattice.sage (updated)
# (c) 2025. Toolkit for ranks and discriminants of Hodge lattices on Fermat hypersurfaces.
# SageMath >= 9.0

from sage.all import *
from itertools import product

############################################################
# Part A: Rank via Shioda–Ran (count |B_m^n|)              #
############################################################

def _units_mod(m):
    return [t for t in range(1, m) if gcd(t, m) == 1]

def _lt_table(m):
    # table[t][r] = < t * (r+1) > in {1,...,m-1}
    tbl = {}
    for t in _units_mod(m):
        row = []
        for r in range(1, m):
            tr = (t*r) % m
            if tr == 0:
                tr = m
            row.append(tr)
        tbl[t] = row
    return tbl

def hodge_rank_via_B(m, n):
    """
    Return 1 + |B_m^n| for even n (n = 2p), where B_m^n is Shioda's set of characters
    that contribute to H^{p,p} in primitive cohomology.
    Implemented by solving the multiplicity-vector ILP:
        sum_r c_r = k = n+2, and sum_r <t*r> c_r = (p+1)m for all t in (Z/mZ)^*.
    For small n (e.g. n=4) this exhaustive enumeration is extremely fast.
    """
    if n % 2 != 0:
        raise ValueError("n must be even (n = 2p).")
    p = n // 2
    k = n + 2
    units = _units_mod(m)
    lt = _lt_table(m)
    target = (p+1) * m

    # Exhaustive enumeration of nonnegative integer tuples (c_1,...,c_{m-1}) summing to k.
    total = Integer(0)
    c = [0] * (m-1)

    def rec(i, remaining):
        nonlocal total
        if i == m-2:
            c[i] = remaining
            # Check linear conditions for all units
            for t in units:
                s = sum(lt[t][r-1] * c[r-1] for r in range(1, m))
                if s != target:
                    return
            # Valid solution
            multin = Integer(factorial(k))
            for r in range(m-1):
                multin //= factorial(c[r])
            total += multin
            return
        for val in range(remaining + 1):
            c[i] = val
            rec(i+1, remaining - val)

    rec(0, k)
    return Integer(1) + total  # add the non-primitive class h^p

############################################################
# Part B1: Linear cycles, Gram matrix, and SNF (fourfolds) #
############################################################

def _perfect_matchings(indices):
    """
    All pairings (perfect matchings) of a finite set of even cardinality.
    For n=4 (|indices|=6) this returns 15 matchings.
    """
    indices = list(indices)
    if not indices:
        return [()]
    i0 = indices[0]
    rest = indices[1:]
    res = []
    for j, ij in enumerate(rest):
        pair = (i0, ij)
        remaining = rest[:j] + rest[j+1:]
        for tail in _perfect_matchings(remaining):
            res.append((pair,) + tail)
    return res

def _linear_cycle_from_matching_and_a(matching, a_tuple, m):
    """
    Represent a linear cycle L as a list of equations of the form x_i - zeta^k x_j = 0.
    Input:
      - matching: tuple of pairs (i,j) pairing coordinates
      - a_tuple:  tuple (a_0,...,a_p) with entries in {0,...,m-1}
      - m: degree (order of root of unity)
    """
    if len(matching) != len(a_tuple):
        raise ValueError("Length mismatch between matching and a_tuple")
    return tuple( (min(matching[t]), max(matching[t]), a_tuple[t] % m) for t in range(len(matching)) )

def _intersection_dim_in_P(cycle1, cycle2, m):
    """
    Compute dim( L1 ∩ L2 ) in projective space P^{n+1} where n=2p, by linear algebra over Q(zeta_m).
    Each cycle is a tuple of triples (i,j,a) encoding equations x_i - zeta^a x_j = 0.
    """
    p = len(cycle1) - 1  # since cycles live in P^{2p+1}
    N = 2*p + 2          # number of homogeneous coordinates = n+2
    K = CyclotomicField(m)
    z = K.gen()
    rows = []
    for (i,j,a) in list(cycle1) + list(cycle2):
        row = [K(0)] * N
        row[i] = K(1)
        row[j] = - z**a  # equation: x_i - z^a x_j = 0
        rows.append(row)
    M = Matrix(K, rows)
    r = M.rank()
    # The vector space ambient is K^N; projectivize => dimension is (N - r - 1)
    return (N - r - 1)

def linear_cycles_intersection_number(cycle1, cycle2, m):
    """
    Return the intersection number (Pi · Pj)_X for X = Fermat fourfold of degree m,
    using the closed formula (AMV Eq.(6)):
         (1 - (-m + 1)^{r+1}) / m,  where r = dim(Pi ∩ Pj) in projective space.
    """
    r = _intersection_dim_in_P(cycle1, cycle2, m)
    return ZZ((1 - Integer(-m + 1)**(r+1)) // m) if r >= -1 else ZZ(0)

def _all_linear_cycles_fourfold(m, max_matchings=None):
    """
    Enumerate linear P^2's on X^4_m. There are 15 pairings and m^3 choices of parameters.
    Use max_matchings to truncate (for experiments).
    Each cycle is encoded as a tuple of triples (i,j,a) with i<j.
    """
    indices = tuple(range(6))  # 0..5
    matchings = _perfect_matchings(indices)
    if max_matchings is not None:
        matchings = matchings[:max_matchings]
    zetas = range(m)  # exponents 0..m-1 for zeta_m^a
    cycles = []
    for mat in matchings:
        for a_tuple in product(zetas, repeat=3):
            cycles.append(_linear_cycle_from_matching_and_a(mat, a_tuple, m))
    return cycles

def _primitive_IP_pair(cycle_i, cycle_j, cycle_base, m):
    """Primitive intersection: (i - base)·(j - base) in X (integer)."""
    Iij = linear_cycles_intersection_number(cycle_i, cycle_j, m)
    Iib = linear_cycles_intersection_number(cycle_i, cycle_base, m)
    Ibj = linear_cycles_intersection_number(cycle_base, cycle_j, m)
    Ibb = linear_cycles_intersection_number(cycle_base, cycle_base, m)
    return Iij - Iib - Ibj + Ibb

def smith_invariants_from_matrix(A):
    """Return invariant factors (nonzero diagonal entries) of the Smith normal form of A."""
    U, S, V = A.smith_form()
    diag = [ZZ(S[i,i]) for i in range(min(S.nrows(), S.ncols()))]
    return [d for d in diag if d != 0]

def adjust_primitive_to_full(invariants, d):
    """
    Adjust invariant factors of the primitive lattice to the full Hodge lattice:
    decrease the multiplicity of 1 by two and append d (cf. AMV, §3).
    """
    inv = list(invariants)
    ones = inv.count(ZZ(1))
    if ones < 2:
        raise ValueError("Not enough 1's to adjust to full lattice; check input invariants.")
    # remove two 1's
    removed = 0
    inv2 = []
    for a in inv:
        if a == 1 and removed < 2:
            removed += 1
        else:
            inv2.append(a)
    inv2.append(ZZ(d))
    inv2.sort()
    return inv2

def linear_cycles_SNF_fourfold(m, max_matchings=None, target_rank=None, seed=0):
    """
    Master routine for fourfolds (n=4):
      1) Build a pool of linear P^2 cycles.
      2) Form the primitive Gram matrix A1 (dense; caution with memory).
      3) Extract a full-rank principal submatrix greedily (target_rank optional).
      4) Compute its Smith normal form invariants.
      5) Return invariants for the primitive lattice, and the adjusted list for the full lattice.
    WARNING: Steps (2–4) can be heavy. Prefer the streaming driver in the separate script.
    """
    set_random_seed(seed)
    cycles = _all_linear_cycles_fourfold(m, max_matchings=max_matchings)
    # Dense build (for small experiments only)
    N = len(cycles)
    base = cycles[0]
    A = Matrix(ZZ, N-1, N-1)
    for i in range(1, N):
        for j in range(1, N):
            A[i-1, j-1] = _primitive_IP_pair(cycles[i], cycles[j], base, m)
    # Estimate target rank if not given
    if target_rank is None:
        target_rank = int(hodge_rank_via_B(m, 4) - 1)  # primitive rank
    # Greedy principal submatrix
    chosen = []
    B = Matrix(QQ, 0, 0)
    for idx in range(N-1):
        cand = A.matrix_from_rows_and_columns(chosen + [idx], chosen + [idx])
        if cand.rank() > B.rank():
            chosen.append(idx)
            B = cand
        if B.rank() >= target_rank:
            break
    prim_inv = smith_invariants_from_matrix(B.change_ring(ZZ))
    full_inv = adjust_primitive_to_full(prim_inv, m)
    disc = prod(full_inv)
    return dict(
        num_cycles=len(cycles),
        chosen=len(chosen),
        primitive_invariants=prim_inv,
        full_invariants=full_inv,
        discriminant=disc
    )
