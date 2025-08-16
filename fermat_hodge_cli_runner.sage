
# fermat_hodge_cli_runner.sage
# Run with:    sage fermat_hodge_cli_runner.sage --m all --b1 --b2-m5 --assert-m5 --log-every 25
#
# This CLI orchestrates:
#   - Rank via Shioda–Ran (Algorithm A)
#   - B1: Linear-cycles route (streaming, greedy basis growth + SNF)
#   - B2: Vanishing-cycles + projectors route (exact primitive Gram + SNF)
#
# Defaults are conservative to avoid very heavy runs:
#   * B2 (projectors) is enabled by default ONLY for m=5 (n=4),
#     because for m=7,9 the primitive Gram is large.
#   * B1 is enabled by default.
#
# Options:
#   --m {5,7,9,all}          which degree(s) to run (default: all)
#   --n 4                    dimension is fixed to 4 in this runner (middle degree), but kept as a flag
#   --b1/--no-b1             enable/disable B1 (default: on)
#   --b2-m5/--no-b2-m5       enable/disable B2 for m=5 (default: on)
#   --b2-others              also run B2 for m=7,9 (warning: heavy)
#   --compute-snf-b1         compute SNF for B1 (default: True for m=5, False otherwise unless set)
#   --max-matchings K        limit number of pairings used in B1 (default: None = use all)
#   --target-rank R          primitive target rank (default: inferred from Shioda–Ran)
#   --log-every k            progress logging interval for B1 (default: 25)
#   --assert-m5              assert published m=5 invariants (exits nonzero on mismatch)
#
from sage.all import *
import sys, time

# Load libraries (must be in the same folder)
try:
    load('fermat_hodge_lattice.sage')         # B1 and rank
    load('fermat_hodge_B2_projectors.sage')   # B2
    load('driver_fourfold_rank_logger.sage')  # streaming B1 growth + optional SNF
except Exception as e:
    print("Error: could not load required .sage files. Make sure they are in the current directory.")
    print(e)
    sys.exit(1)

def known_m5_full_invariants():
    # AMV Table 1 for (n,d)=(4,5): 1^166 · 5^174 · 25^54 · 125^7
    return [ZZ(1)]*166 + [ZZ(5)]*174 + [ZZ(25)]*54 + [ZZ(125)]*7

def summarize_invariants(inv):
    # Compact summary "a^k * b^l * ..."
    from collections import Counter
    c = Counter(inv)
    terms = []
    for val in sorted(c.keys()):
        exp = c[val]
        terms.append(f"{val}^{exp}" if exp != 1 else f"{val}")
    return " * ".join(terms)

def run_case(m=5, n=4, use_b1=True, use_b2_m5=True, use_b2_others=False,
             compute_snf_b1=None, max_matchings=None, target_rank=None,
             log_every=25, assert_m5=False, seed=1):
    print("="*80)
    print(f"Fermat fourfold: degree m={m}, dimension n={n} (middle cohomology in H^4)")
    print("="*80)

    # Rank via Shioda–Ran
    R = hodge_rank_via_B(m, n)
    primitive_target = R - 1
    print(f"[Rank/A] Shioda–Ran gives full rank = {R}, primitive target = {primitive_target}.")

    # Decide which routes to run
    run_b1 = use_b1
    run_b2 = (use_b2_m5 if m == 5 else use_b2_others)

    # B2: vanishing cycles + projectors (primitive Gram + SNF)
    prim_inv_B2 = full_inv_B2 = disc_B2 = None
    if run_b2:
        t0 = time.time()
        print(f"[B2] Building primitive Gram via vanishing cycles/projectors (this may take a while)...")
        G = primitive_gram_vanishing(m, n)
        print(f"[B2] Primitive Gram: size = {G.nrows()} x {G.ncols()} (expected {primitive_target}).")
        if G.nrows() != primitive_target:
            print("[B2] WARNING: primitive Gram size did not match target; proceeding anyway.")
        print(f"[B2] Computing SNF invariants...")
        prim_inv_B2 = smith_invariants_from_matrix(G)
        full_inv_B2  = full_invariants_from_primitive(prim_inv_B2, m)
        disc_B2 = prod(full_inv_B2)
        dt = time.time() - t0
        print(f"[B2] Done in {dt:.1f}s. Full-lattice inv summary: {summarize_invariants(full_inv_B2)}")
        print(f"[B2] Discriminant (full lattice): {disc_B2}")

    # B1: linear cycles (streaming greedy + optional SNF)
    res_b1 = None
    if run_b1:
        if compute_snf_b1 is None:
            compute_snf_b1 = (m == 5)
        print(f"[B1] Streaming greedy growth (log every {log_every}). "
              f"compute_snf={compute_snf_b1}, max_matchings={max_matchings}, target_rank={target_rank}")
        res_b1 = grow_basis_and_log(m, max_matchings=max_matchings,
                                    target_rank=(primitive_target if target_rank is None else target_rank),
                                    log_every=log_every, shuffle=True,
                                    compute_snf=compute_snf_b1, seed=seed)
        print("[B1] Summary:", res_b1)
        if compute_snf_b1:
            print(f"[B1] Full-lattice inv summary: {summarize_invariants(res_b1['full_invariants'])}")
            print(f"[B1] Discriminant (full lattice): {res_b1['discriminant']}")

    # Compare B1 and B2 if both available
    if run_b1 and run_b2 and res_b1 is not None and 'full_invariants' in res_b1:
        same_prim = (res_b1.get('full_invariants') is not None and prim_inv_B2 is not None and 
                     res_b1['primitive_rank'] == len(prim_inv_B2) )
        # Compare lists exactly
        same_full = (res_b1.get('full_invariants') == full_inv_B2)
        print(f"[COMPARE] primitive lengths equal? {same_prim} (B1 rank {res_b1.get('primitive_rank')}, B2 size {len(prim_inv_B2) if prim_inv_B2 else 'N/A'})")
        print(f"[COMPARE] full invariant lists equal? {same_full}")
        if not same_full:
            print("[COMPARE] B1 full invariants (summary):", summarize_invariants(res_b1['full_invariants']))
            print("[COMPARE] B2 full invariants (summary):", summarize_invariants(full_inv_B2))

    # Assert published invariants for m=5 if requested
    if assert_m5 and m == 5:
        if res_b1 is None or 'full_invariants' not in res_b1:
            print("[ASSERT m=5] ERROR: B1 full invariants not available; cannot assert.")
            sys.exit(2)
        expected = known_m5_full_invariants()
        got = list(res_b1['full_invariants'])
        ok = (got == expected)
        print(f"[ASSERT m=5] match with published list? {ok}")
        if not ok:
            print("[ASSERT m=5] Expected summary:", summarize_invariants(expected))
            print("[ASSERT m=5] Got      summary:", summarize_invariants(got))
            sys.exit(3)

def parse_args(argv):
    # Minimal positional/flag parsing compatible with Sage
    # Use a very light parser to avoid argparse edge cases in Sage mode
    opts = {
        'm_list': [5,7,9],
        'n': 4,
        'b1': True,
        'b2_m5': True,
        'b2_others': False,
        'compute_snf_b1': None,
        'max_matchings': None,
        'target_rank': None,
        'log_every': 25,
        'assert_m5': False,
        'seed': 1,
    }
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ('--m', '-m'):
            i += 1; val = argv[i]
            if val == 'all':
                opts['m_list'] = [5,7,9]
            else:
                opts['m_list'] = [int(x) for x in val.split(',')]
        elif a == '--n':
            i += 1; opts['n'] = int(argv[i])
        elif a == '--b1':
            opts['b1'] = True
        elif a == '--no-b1':
            opts['b1'] = False
        elif a == '--b2-m5':
            opts['b2_m5'] = True
        elif a == '--no-b2-m5':
            opts['b2_m5'] = False
        elif a == '--b2-others':
            opts['b2_others'] = True
        elif a == '--compute-snf-b1':
            i += 1; val = argv[i].lower()
            if val in ('true','1','yes','y'):  opts['compute_snf_b1'] = True
            if val in ('false','0','no','n'):  opts['compute_snf_b1'] = False
        elif a == '--max-matchings':
            i += 1; opts['max_matchings'] = int(argv[i])
        elif a == '--target-rank':
            i += 1; opts['target_rank'] = int(argv[i])
        elif a == '--log-every':
            i += 1; opts['log_every'] = int(argv[i])
        elif a == '--assert-m5':
            opts['assert_m5'] = True
        elif a == '--seed':
            i += 1; opts['seed'] = int(argv[i])
        i += 1
    return opts

if __name__ == "__main__":
    # Drop the script name
    argv = sys.argv[1:]
    opts = parse_args(argv)

    for m in opts['m_list']:
        run_case(m=m, n=opts['n'],
                 use_b1=opts['b1'],
                 use_b2_m5=opts['b2_m5'],
                 use_b2_others=opts['b2_others'],
                 compute_snf_b1=opts['compute_snf_b1'],
                 max_matchings=opts['max_matchings'],
                 target_rank=opts['target_rank'],
                 log_every=opts['log_every'],
                 assert_m5=opts['assert_m5'],
                 seed=opts['seed'])
