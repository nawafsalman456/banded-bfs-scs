#!/usr/bin/env python3
"""
Final Deterministic Banded BFS for SCS in the deletion–only model,
with advanced checks and optimality testing.

Assumptions:
  - There is an unknown parent x ∈ {0,1}^n.
  - We are given k binary traces, each obtained from x by at most t deletions.
  - Our goal is to compute the shortest common supersequence (SCS) z such that each trace is a subsequence of z.
  
We use a banded BFS over states defined by a tuple of pointers (one per trace). 
A move outputs one letter. For each trace, if the letter equals its next required symbol, we advance its pointer;
otherwise, we “skip” that trace (counted as a deletion). We enforce that for every active trace j,
   (output_length – pointer_j) ≤ t.
We also restrict the state space by requiring that the difference between the maximum and minimum pointers is at most L_band = L_factor * t.

This module includes:
  1. calc_scs_bfs: the banded BFS algorithm.
  2. calc_scs_dp: a DP SCS finder.
  3. A suite of tests (basic and advanced) that check correctness, optimality, and timing.

Authors: Nawaf Salman, Nahel Awidat
"""

import random
import time
import time
import matplotlib.pyplot as plt
from collections import deque
from itertools import permutations, combinations

DEBUG_MODE = False  # change to True to run in debug mode (enables more prints)

def debug_print(*args, **kwargs):
    if DEBUG_MODE:
        print(*args, **kwargs)

# -------------------------------------------------
# Helper: Check if 'small' is a subsequence of 'large'
# -------------------------------------------------
def is_subsequence(small, large):
    it = iter(large)
    return all(c in it for c in small)

# -------------------------------------------------
# Our Banded BFS Algorithm for SCS
# -------------------------------------------------
def calc_scs_bfs(sequences, t):
    """
    Compute the shortest common supersequence (SCS) of the input binary sequences
    under the deletion-only model. Each sequence in `sequences` is a string of '0' and '1'.
    
    Parameter:
        sequences: list of binary strings (each obtained from the same x with at most t deletions)
        t: maximum allowed number of deletions for each sequence.
        
    Returns:
        (scs, cost) where scs is the SCS (a string) and cost is its length.
    """
    k = len(sequences)
    # Precompute lengths.
    lengths = [len(seq) for seq in sequences]
    
    # A state is an k-tuple of indices (0-indexed). Start state is (0,0,...,0).
    start = tuple(0 for _ in range(k))
    goal = tuple(lengths[i] for i in range(k))
    
    # Each state has an associated cost (length of SCS built so far), parent, and char that was appended.
    # We'll use a dictionary: state -> (cost, parent_state, char)
    dp = {start: (0, None, None)}
    q = deque([start])
    
    # Deletion invariant: for every active sequence i (i.e. p_i < len(seq)), cost - p_i <= t.
    def valid_state(state, cost):
        for i in range(k):
            if state[i] < lengths[i]:
                if cost - state[i] > t:
                    return False
        return True

    # Reconstruct the supersequence from dp.
    def reconstruct(state):
        chars = []
        while dp[state][1] is not None:
            cost, parent, char = dp[state]
            chars.append(char)
            state = parent
        return "".join(reversed(chars))
    
    # BFS: each transition increases cost by 1.
    while q:
        state = q.popleft()
        cost, _, _ = dp[state]
        
        # Check if goal reached.
        if state == goal:
            scs = reconstruct(state)
            return reconstruct(state), cost
        
        # Determine candidate letters from active sequences.
        # R is the set of letters that appear as the next character in at least one active sequence.
        R = set()
        for i in range(k):
            if state[i] < lengths[i]:
                R.add(sequences[i][state[i]])
        # If R is empty (should not happen because if all are exhausted, state==goal), skip.
        if not R:
            continue
        
        # For each candidate letter, compute the new state.
        for letter in R:
            new_state = []
            for i in range(k):
                # If sequence i is active and its next letter equals 'letter', advance pointer.
                if state[i] < lengths[i] and sequences[i][state[i]] == letter:
                    new_state.append(state[i] + 1)
                else:
                    new_state.append(state[i])
            new_state = tuple(new_state)
            new_cost = cost + 1
            
            # Check deletion invariant: for each active sequence, new_cost - new_state[i] <= t.
            if not valid_state(new_state, new_cost):
                continue
            
            # If new_state has not been seen or we found a shorter way, record it.
            if new_state not in dp or dp[new_state][0] > new_cost:
                dp[new_state] = (new_cost, state, letter)
                q.append(new_state)
                
    # If goal is unreachable (should not happen under our assumptions), return None.
    return None, float('inf')

# -------------------------------------------------
# Known DP SCS
# -------------------------------------------------

def scs_two(X, Y):
    """
    Returns the Shortest Common Supersequence (SCS) of two words using a DP table 
    + backtracking reconstruction. If one string is a subsequence of the other, 
    it returns the longer one immediately.
    """
    if is_subsequence(X, Y):
        return Y
    if is_subsequence(Y, X):
        return X

    m, n = len(X), len(Y)
    # DP table for lengths of SCS
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if X[i - 1] == Y[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = min(dp[i - 1][j], dp[i][j - 1]) + 1

    # Reconstruct SCS from dp table
    i, j = m, n
    result = []
    while i > 0 and j > 0:
        if X[i - 1] == Y[j - 1]:
            result.append(X[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] < dp[i][j - 1]:
            result.append(X[i - 1])
            i -= 1
        else:
            result.append(Y[j - 1])
            j -= 1
    # Append leftovers
    while i > 0:
        result.append(X[i - 1])
        i -= 1
    while j > 0:
        result.append(Y[j - 1])
        j -= 1

    return "".join(reversed(result))

def verify_scs(scs_result, words):
    """Verifies that all words are subsequences of the computed SCS."""
    return all(is_subsequence(word, scs_result) for word in words)

def calc_scs_dp(seqs):
    """
    Find the globally shortest SCS for a list of words by:
      1. Checking if any single word is already a supersequence of all.
      2. Otherwise, for each permutation of seqs, merge from left to right
         using 'scs_two' and pick the overall shortest.
    """
    words = list(seqs)
    # Quick check: if any word is already a supersequence of all, pick the shortest such word.
    superseq_candidates = [w for w in words if all(is_subsequence(s, w) for s in words)]
    if superseq_candidates:
        return min(superseq_candidates, key=len)

    best_scs = None
    # Try all permutations of the input words:
    for perm in permutations(words):
        merged = perm[0]
        for i in range(1, len(perm)):
            merged = scs_two(merged, perm[i])
        # Keep track of the globally shortest
        if best_scs is None or len(merged) < len(best_scs):
            best_scs = merged

    # Final check
    if not verify_scs(best_scs, seqs):
        raise ValueError("calc_scs_dp failed to produce a valid SCS.")

    return best_scs

# -------------------------------------------------
# Test Data Generation Utilities
# -------------------------------------------------
def generate_random_x(n):
    """Generate a random binary string of length n."""
    return "".join(random.choice("01") for _ in range(n))

def generate_trace(x, t, delete_exactly=True):
    """
    Generate a trace from x by deleting up to t symbols.
    If delete_exactly is True, exactly t symbols are deleted; otherwise, a random number in [0,t] is deleted.
    Returns (trace, deletion_set).
    """
    n = len(x)
    indices = list(range(n))
    if delete_exactly:
        num_del = t
    else:
        num_del = random.randint(0, t)
    dset = set(random.sample(indices, num_del))
    trace = "".join(x[i] for i in range(n) if i not in dset)
    return trace, dset

def generate_traces(x, k, t, delete_exactly=True, ensure_nonuniversal=True):
    """
    Generate k traces from x by deleting up to t symbols (exactly t if delete_exactly=True).
    If ensure_nonuniversal is True, ensure that for every position in x at least one trace keeps that symbol.
    Returns a list of traces.
    """
    n = len(x)
    traces = []
    deletion_sets = []
    for _ in range(k):
        trace, dset = generate_trace(x, t, delete_exactly=delete_exactly)
        traces.append(trace)
        deletion_sets.append(dset)
    if ensure_nonuniversal:
        for pos in range(n):
            if all(pos in ds for ds in deletion_sets):
                j = random.randrange(k)
                deletion_sets[j].remove(pos)
                traces[j] = "".join(x[i] for i in range(n) if i not in deletion_sets[j])
    return traces


# -------------------------------------------------
# Basic Checks
# -------------------------------------------------
def test_scs_basic(seqs, t, expected=None):
    scs, cost = calc_scs_bfs(seqs, t)
    debug_print("Basic Test: sequences =", seqs)
    debug_print("  Computed SCS:", scs, "with length", cost)
    for seq in seqs:
        if not is_subsequence(seq, scs):
            print("  ERROR: {} is not a subsequence of {}".format(seq, scs))
            exit(-1)
    if expected is not None and scs != expected:
        print("  Note: expected SCS =", expected)
        exit(-1)
    return True

# -------------------------------------------------
# Advanced Checks
# -------------------------------------------------
def advanced_checks():
    debug_print("=== Advanced Checks ===")
    
    # Test 1: Small instance optimality test.
    debug_print("Test 1: small instance optimality check")
    x_small = generate_random_x(20)
    r_small = 3
    t_small = 2
    traces1 = generate_traces(x_small, r_small, t_small, delete_exactly=True, ensure_nonuniversal=True)
    debug_print("Parent x_small:", x_small)
    debug_print("Traces:", traces1)
    # start = time.time()
    scs1, cost1 = calc_scs_bfs(traces1, t_small)
    # elapsed = time.time() - start
    # print("Elapsed time: {:.4f} seconds".format(elapsed))
    scs_dp = calc_scs_dp(traces1)
    scs_dp_len = len(scs_dp)
    debug_print("Banded BFS SCS length:", cost1)
    debug_print("Known DP SCS length:", scs_dp_len)
    if cost1 == scs_dp_len:
        # lenght comparision is good enough:
        #   - scs_dp function returns a valid common supersequence (else it will fail before returning)
        #   - if scs_dp_len == cost1, it means calc_scs_bfs returned a valid SCS
        # we can't test the SCS directly because 
        debug_print("Test 1 PASSED: Optimality confirmed.\n")
    else:
        print("Test 1 FAILED: SCS lenght mismatch.\n")
        print("traces   - ", traces1)
        print("expected - ", scs_dp)
        print("actual   - ", scs1)
        exit(-1)
    
    # Test 2: Varying deletion counts (some sequences delete less than t).
    debug_print("Test 2: Varying deletion counts (delete up to t).")
    x2 = generate_random_x(50)
    r2 = 3
    t2 = 5
    traces2 = []
    for _ in range(r2):
        trace, _ = generate_trace(x2, t2, delete_exactly=False)
        traces2.append(trace)
    scs2, cost2 = calc_scs_bfs(traces2, t2)
    debug_print("Parent x2:", x2)
    debug_print("Traces:", traces2)
    debug_print("Computed SCS:", scs2, "length:", cost2)
    if verify_scs(scs2, traces2):
        debug_print("Test 2 PASSED: All traces are subsequences.\n")
    else:
        print("Test 2 FAILED.\n")
        exit(-1)
    
    # Test 3: Universal deletion (SCS not equal to x).
    debug_print("Test 3: Universal deletion (forcing SCS != x).")
    x3 = generate_random_x(100)
    r3 = 3
    t3 = 10
    traces3 = []
    # Force a particular position (pos0) to be deleted in every trace.
    pos0 = random.randint(0,100)
    for _ in range(r3):
        trace, dset = generate_trace(x3, t3, delete_exactly=True)
        dset.add(pos0)
        new_trace = "".join(x3[i] for i in range(len(x3)) if i not in dset)
        traces3.append(new_trace)
    scs3, cost3 = calc_scs_bfs(traces3, t3)
    debug_print("Parent x3:", x3)
    debug_print("Traces:", traces3)
    debug_print("Computed SCS:", scs3, "length:", cost3)
    if verify_scs(scs3, traces3):
        debug_print("Test 3 PASSED: All traces are subsequences.\n")
    else:
        print("Test 3 FAILED.\n")
        exit(-1)
    
    # Test 4: Timing test on moderately long instance.
    debug_print("Test 4: Timing test on long instance.")
    n4 = 4000
    x4 = generate_random_x(n4)
    r4 = 3
    t4 = 5
    traces4 = generate_traces(x4, r4, t4, delete_exactly=True, ensure_nonuniversal=True)
    debug_print("Start BFS scs search...")
    start = time.time()
    scs4, cost4 = calc_scs_bfs(traces4, t4)
    elapsed = time.time() - start
    debug_print("Elapsed time: {:.4f} seconds".format(elapsed))
    
    debug_print("Start DP scs search...")
    start = time.time()
    _ = calc_scs_dp(traces4)
    elapsed = time.time() - start
    debug_print("Elapsed time: {:.4f} seconds".format(elapsed))

    if verify_scs(scs4, traces4):
        debug_print("Test 4 PASSED: All traces are subsequences.\n")
    else:
        print("Test 4 FAILED.\n")
        exit(-1)
    

def experiment_timing():
    
    # 1) n in [10..1000], k=3, t=5 => measure BFS time and DP time, plot times vs n.
    Nvals = list(range(10,1000,10))
    BFS_times_n = []
    DP_times_n = []
    k_fixed = 3
    t_fixed = 5
    for n in Nvals:
        debug_print(f"\nExperiment: n={n}, k={k_fixed}, t={t_fixed}")
        x = generate_random_x(n)
        # Each trace => length n-t
        traces = generate_traces(x, k_fixed, t_fixed)
        
        # BFS timing
        start = time.time()
        scs_bfs_res, bfs_len = calc_scs_bfs(traces, t_fixed)
        bfs_time = time.time() - start
        
        # DP timing
        start = time.time()
        scs_dp_res = calc_scs_dp(traces)
        dp_time = time.time() - start
        
        BFS_times_n.append(bfs_time)
        DP_times_n.append(dp_time)
        
        # Optional correctness check:
        if not verify_scs(scs_bfs_res, traces):
            print(f"ERROR: BFS SCS is not valid for n={n}!")
            exit(-1)
            # can break or continue
        if not verify_scs(scs_dp_res, traces):
            print(f"ERROR: DP SCS is not valid for n={n}!")
            exit(-1)
        debug_print(f" Banded BFS time={bfs_time:.4f}, length={bfs_len}  DP time={dp_time:.4f}, length={len(scs_dp_res)}")
    
    plt.figure(figsize=(8,6))
    plt.plot(Nvals, BFS_times_n, label="Banded BFS", marker='o')
    plt.plot(Nvals, DP_times_n, label="DP merges", marker='s')
    plt.xlabel("n")
    plt.ylabel("Time (seconds)")
    plt.title(f"Timing vs. n (k={k_fixed}, t={t_fixed})")
    plt.legend()
    plt.savefig("scs_calc_times_vs_n")
    
    # 2) k in [3..10], fix n=20, t=5 => measure BFS time and DP time, plot times vs k.
    kvals = list(range(3,11))
    BFS_times_k = []
    DP_times_k = []
    n_fixed = 20
    t_fixed = 3
    
    x = generate_random_x(n_fixed)
    for k in kvals:
        debug_print(f"\nExperiment: n={n_fixed}, k={k}, t={t_fixed}")
        traces = generate_traces(x, k, t_fixed)
        
        start = time.time()
        scs_bfs_res, bfs_len = calc_scs_bfs(traces, t_fixed)
        bfs_time = time.time() - start
        
        start = time.time()
        scs_dp_res = calc_scs_dp(traces)
        dp_time = time.time() - start
        
        BFS_times_k.append(bfs_time)
        DP_times_k.append(dp_time)
        
        if not verify_scs(scs_bfs_res, traces):
            print(f"ERROR: BFS SCS is not valid for k={k}!")
            exit(-1)
        if not verify_scs(scs_dp_res, traces):
            print(f"ERROR: DP SCS is not valid for k={k}!")
            exit(-1)
        debug_print(f" Banded BFS time={bfs_time:.4f}, length={bfs_len}  DP time={dp_time:.4f}, length={len(scs_dp_res)}")
    
    plt.figure(figsize=(8,6))
    plt.plot(kvals, BFS_times_k, label="Banded BFS", marker='o')
    plt.plot(kvals, DP_times_k, label="DP merges", marker='s')
    plt.xlabel("k")
    plt.ylabel("Time (seconds)")
    plt.title(f"Timing vs. k (n={n_fixed}, t={t_fixed})")
    plt.legend()
    plt.savefig("scs_calc_times_vs_k")
    

# -------------------------------------------------
# Main Testing Routine
# -------------------------------------------------
if __name__ == "__main__":
    # Run a few basic tests
    print("start basic tests ...")
    test_scs_basic(["10101", "1001"], t=1, expected="10101")
    test_scs_basic(["1100", "1100"], t=0, expected="1100")
    test_scs_basic(["1010", "1000", "1100"], t=1)
    print("Basic tests PASSED.")
    
    # Run more advanced checks
    print("start advanced tests ...")
    advanced_checks()
    print("Advanced tests PASSED.")
    
    # Run timing experiments
    print("start timing experiments (this will take few minutes because of the non-optimal DP algorithm) ...")
    experiment_timing()
    print("Timing experiments DONE.")
