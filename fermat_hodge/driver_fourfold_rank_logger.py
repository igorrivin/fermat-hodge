
# driver_fourfold_rank_logger.sage
# A streaming driver to watch the primitive rank grow for Fermat fourfolds (m=7,9),
# and (optionally) compute the Smith normal form at the end.
#
# Place this file in the same folder as 'fermat_hodge_lattice.sage' and run in Sage:
#     sage -python driver_fourfold_rank_logger.sage
# or within a Sage REPL:
#     load('driver_fourfold_rank_logger.sage')
#
from sage.all import *
import time

# Load the library
load('fermat_hodge_lattice.sage')

def grow_basis_and_log(m, max_matchings=None, target_rank=None, log_every=1, shuffle=False, compute_snf=False, seed=0):
    """
    Streaming greedy construction of a primitive Gram matrix for linear P^2's on X_m^4.
    Logs rank increases so you can 'watch' it tick up to the target.
    Returns a dict with timing, chosen count, and (optionally) SNF info.
    """
    if target_rank is None:
        target_rank = int(hodge_rank_via_B(m, 4) - 1)
    if log_every < 1:
        log_every = 1

    cycles = _all_linear_cycles_fourfold(m, max_matchings=max_matchings)
    if shuffle:
        set_random_seed(seed)
        from random import shuffle as pyshuffle
        pyshuffle(cycles)

    base = cycles[0]
    chosen = []       # indices into 'cycles' (excluding base)
    B = Matrix(QQ, 0, 0)  # current primitive Gram (r x r)

    t0 = time.time()
    last_logged = 0
    added = 0

    # Precompute base self intersection (constant)
    Ibb = linear_cycles_intersection_number(base, base, m)

    def primitive_ip_with_chosen(idx):
        """Return (w, s) where w_j = (Pi - base)·(Pj - base), s = self primitive for Pi."""
        Pi = cycles[idx]
        # vector against chosen
        w = []
        for jdx in chosen:
            Pj = cycles[jdx]
            val = _primitive_IP_pair(Pi, Pj, base, m)
            w.append(ZZ(val))
        # self term
        s = _primitive_IP_pair(Pi, Pi, base, m)
        return vector(QQ, w), QQ(s)

    # Iterate candidates
    for idx in range(1, len(cycles)):
        w, s = primitive_ip_with_chosen(idx)
        # Form trial augmented matrix
        if B.nrows() == 0:
            trial_rank = 1 if s != 0 else 0
            rank_increased = (trial_rank > 0)
        else:
            # augment B with new row/col
            top = B.augment(matrix(QQ, B.nrows(), 1, w))
            bottom = matrix(QQ, 1, B.ncols(), w).augment(matrix(QQ, 1, 1, [s]))
            B_trial = top.stack(bottom)
            rank_increased = (B_trial.rank() > B.rank())

        if rank_increased:
            # Commit the augmentation
            if B.nrows() == 0:
                B = matrix(QQ, 1, 1, [s])
            else:
                top = B.augment(matrix(QQ, B.nrows(), 1, w))
                bottom = matrix(QQ, 1, B.ncols(), w).augment(matrix(QQ, 1, 1, [s]))
                B = top.stack(bottom)
            chosen.append(idx)
            added += 1
            if added % log_every == 0 or added <= 5 or added == target_rank:
                elapsed = time.time() - t0
                print(f"[m={m}] primitive rank {B.rank()}/{target_rank}  (+1)  after {elapsed:8.1f}s  (added #{added}, candidate idx {idx})")

            if B.rank() >= target_rank:
                break

    out = {
        'm': m,
        'num_cycles_seen': idx+1,
        'chosen': len(chosen),
        'primitive_rank': int(B.rank()),
        'target_rank': target_rank,
        'elapsed_sec': time.time() - t0,
    }

    if compute_snf:
        print(f"[m={m}] Computing Smith normal form on a {B.nrows()}x{B.ncols()} integer matrix...")
        # Convert to ZZ and compute SNF
        Bz = B.change_ring(ZZ)
        prim_inv = smith_invariants_from_matrix(Bz)
        full_inv = adjust_primitive_to_full(prim_inv, m)
        disc = prod(full_inv)
        out.update({
            'primitive_invariants': prim_inv,
            'full_invariants': full_inv,
            'discriminant': disc
        })
        print(f"[m={m}] Done. Discriminant (full lattice) = {disc}")

    return out

def run_demo(m, log_every=10, compute_snf=False, seed=0):
    print(f"=== Demo for Fermat fourfold of degree m={m} ===")
    R = hodge_rank_via_B(m, 4)
    print(f"Target full rank (including h^2) = {R}; primitive target = {R-1}.")
    res = grow_basis_and_log(m, log_every=log_every, compute_snf=compute_snf, seed=seed)
    print("Summary:", res)
    return res

if __name__ == "__main__":
    # Gentle defaults: log every 25 new basis elements; don't compute SNF by default.
    run_demo(7, log_every=25, compute_snf=False, seed=1)
    # Uncomment the next line if you want to try m=9 (can be longer).
    # run_demo(9, log_every=25, compute_snf=False, seed=1)
