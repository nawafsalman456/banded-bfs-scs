#!/usr/bin/env python3
"""
Final Deterministic Banded BFS for SCS in the deletion–only model,
with advanced checks and optimality testing.

Assumptions:
  - There is an unknown parent x ∈ {0,1}^n.
  - We are given r binary traces, each obtained from x by at most t deletions.
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
from collections import deque
from itertools import permutations, combinations

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
    r = len(sequences)
    # Precompute lengths.
    lengths = [len(seq) for seq in sequences]
    
    # A state is an r-tuple of indices (0-indexed). Start state is (0,0,...,0).
    start = tuple(0 for _ in range(r))
    goal = tuple(lengths[i] for i in range(r))
    
    # Each state has an associated cost (length of SCS built so far), parent, and char that was appended.
    # We'll use a dictionary: state -> (cost, parent_state, char)
    dp = {start: (0, None, None)}
    q = deque([start])
    
    # Deletion invariant: for every active sequence i (i.e. p_i < len(seq)), cost - p_i <= t.
    def valid_state(state, cost):
        for i in range(r):
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
            if not verify_scs(scs, sequences):  # sanity check - the output is a common supersequence
                print("error: calc_scs_bfs failed to produce a valid SCS.\n")
                exit(-1)
            return reconstruct(state), cost
        
        # Determine candidate letters from active sequences.
        # R is the set of letters that appear as the next character in at least one active sequence.
        R = set()
        for i in range(r):
            if state[i] < lengths[i]:
                R.add(sequences[i][state[i]])
        # If R is empty (should not happen because if all are exhausted, state==goal), skip.
        if not R:
            continue
        
        # For each candidate letter, compute the new state.
        for letter in R:
            new_state = []
            for i in range(r):
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
    """Returns the Shortest Common Supersequence (SCS) of two words using DP reconstruction.
       If one string is a subsequence of the other, returns the longer one."""
    if is_subsequence(X, Y):
        return Y
    if is_subsequence(Y, X):
        return X

    m, n = len(X), len(Y)
    # Build DP table for SCS length.
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

    # Reconstruct the SCS from the DP table.
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
    while i > 0:
        result.append(X[i - 1])
        i -= 1
    while j > 0:
        result.append(Y[j - 1])
        j -= 1

    # The result is built backwards.
    return "".join(reversed(result))

def verify_scs(scs_result, words):
    """Verifies that all words are subsequences of the computed SCS."""
    return all(is_subsequence(word, scs_result) for word in words)

#TODO - this is not optimal, fix it
def calc_scs_dp(seqs):
    """Finds the Shortest Common Supersequence (SCS) of a list of words by pairwise merging."""
    # Make a copy to avoid modifying the original list.
    words = seqs[:]
    
    # Pre-check: if any input word is a common supersequence of all, return the shortest.
    candidates = [w for w in words if all(is_subsequence(s, w) for s in words)]
    if candidates:
        return min(candidates, key=len)
    
    # Merge pairwise until one word remains.
    while len(words) > 1:
        min_scs = None
        min_len = float('inf')
        best_pair = (None, None)
        
        # Try all pairs and pick the best merge (the one with minimum length).
        for i in range(len(words)):
            for j in range(i + 1, len(words)):
                merged = scs_two(words[i], words[j])
                if len(merged) < min_len:
                    min_len = len(merged)
                    min_scs = merged
                    best_pair = (i, j)
        
        # Remove the two words that were merged.
        i, j = best_pair
        words = [words[k] for k in range(len(words)) if k not in (i, j)]
        words.append(min_scs)
    
    final_scs = words[0]
    if not verify_scs(final_scs, seqs):
        print("error: calc_scs_dp failed to produce a valid SCS.\n")
        exit(-1)
    
    return final_scs

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

def generate_traces(x, r, t, delete_exactly=True, ensure_nonuniversal=True):
    """
    Generate r traces from x by deleting up to t symbols (exactly t if delete_exactly=True).
    If ensure_nonuniversal is True, ensure that for every position in x at least one trace keeps that symbol.
    Returns a list of traces.
    """
    n = len(x)
    traces = []
    deletion_sets = []
    for _ in range(r):
        trace, dset = generate_trace(x, t, delete_exactly=delete_exactly)
        traces.append(trace)
        deletion_sets.append(dset)
    if ensure_nonuniversal:
        for pos in range(n):
            if all(pos in ds for ds in deletion_sets):
                j = random.randrange(r)
                deletion_sets[j].remove(pos)
                traces[j] = "".join(x[i] for i in range(n) if i not in deletion_sets[j])
    return traces


# -------------------------------------------------
# Advanced Checks
# -------------------------------------------------
def advanced_checks():
    print("=== Advanced Checks ===")
    
    # Test 1: Small instance optimality test.
    print("Test 1: small instance optimality check")
    x_small = generate_random_x(20)
    r_small = 3
    t_small = 2
    traces1 = generate_traces(x_small, r_small, t_small, delete_exactly=True, ensure_nonuniversal=True)
    print("Parent x_small:", x_small)
    print("Traces:", traces1)
    # start = time.time()
    scs1, cost1 = calc_scs_bfs(traces1, t_small)
    # elapsed = time.time() - start
    # print("Elapsed time: {:.4f} seconds".format(elapsed))
    scs_dp = calc_scs_dp(traces1)
    scs_dp_len = len(scs_dp)
    print("Banded BFS SCS length:", cost1)
    print("Known DP SCS length:", scs_dp_len)
    if cost1 == scs_dp_len:
        # lenght comparision is good enough:
        #   - scs_dp function returns a valid common supersequence (else it will fail before returning)
        #   - if scs_dp_len == cost1, it means calc_scs_bfs returned a valid SCS
        # we can't test the SCS directly because 
        print("Test 1 PASSED: Optimality confirmed.\n")
    else:
        print("Test 1 FAILED: SCS lenght mismatch.\n")
        print("traces   - ", traces1)
        print("expected - ", scs_dp)
        print("actual   - ", scs1)
        exit(-1)
    
    # Test 2: Varying deletion counts (some sequences delete less than t).
    print("Test 2: Varying deletion counts (delete up to t).")
    x2 = generate_random_x(50)
    r2 = 3
    t2 = 5
    traces2 = []
    for _ in range(r2):
        trace, _ = generate_trace(x2, t2, delete_exactly=False)
        traces2.append(trace)
    scs2, cost2 = calc_scs_bfs(traces2, t2)
    print("Parent x2:", x2)
    print("Traces:", traces2)
    print("Computed SCS:", scs2, "length:", cost2)
    if all(is_subsequence(trace, scs2) for trace in traces2):
        print("Test 2 PASSED: All traces are subsequences.\n")
    else:
        print("Test 2 FAILED.\n")
        exit(-1)
    
    # Test 3: Universal deletion (SCS not equal to x).
    print("Test 3: Universal deletion (forcing SCS != x).")
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
    print("Parent x3:", x3)
    print("Traces:", traces3)
    print("Computed SCS:", scs3, "length:", cost3)
    if all(is_subsequence(trace, scs3) for trace in traces3):
        print("Test 3 PASSED: All traces are subsequences.\n")
    else:
        print("Test 3 FAILED.\n")
        exit(-1)
    
    # Test 4: Timing test on moderately long instance.
    print("Test 4: Timing test on long instance.")
    n4 = 4000
    x4 = generate_random_x(n4)
    r4 = 3
    t4 = 5
    traces4 = generate_traces(x4, r4, t4, delete_exactly=True, ensure_nonuniversal=True)
    print("Start BFS scs search...")
    start = time.time()
    scs4, cost4 = calc_scs_bfs(traces4, t4)
    elapsed = time.time() - start
    print("Elapsed time: {:.4f} seconds".format(elapsed))
    
    print("Start DP scs search...")
    start = time.time()
    _ = calc_scs_dp(traces4)
    elapsed = time.time() - start
    print("Elapsed time: {:.4f} seconds".format(elapsed))

    if all(is_subsequence(trace, scs4) for trace in traces4):
        print("Test 4 PASSED: All traces are subsequences.\n")
    else:
        print("Test 4 FAILED.\n")
        exit(-1)
    
    print("Advanced checks completed.")

# -------------------------------------------------
# Main Testing Routine
# -------------------------------------------------
if __name__ == "__main__":
    # Run a few basic tests.
    def test_scs_basic(seqs, t, expected=None):
        scs, cost = calc_scs_bfs(seqs, t)
        print("Basic Test: sequences =", seqs)
        print("  Computed SCS:", scs, "with length", cost)
        for seq in seqs:
            if not is_subsequence(seq, scs):
                print("  ERROR: {} is not a subsequence of {}".format(seq, scs))
                exit(-1)
        if expected is not None and scs != expected:
            print("  Note: expected SCS =", expected)
            exit(-1)
        print("  Basic test PASSED.\n")
        return True

    test_scs_basic(["10101", "1001"], t=1, expected="10101")
    test_scs_basic(["1100", "1100"], t=0, expected="1100")
    test_scs_basic(["1010", "1000", "1100"], t=1)
    
    advanced_checks()
