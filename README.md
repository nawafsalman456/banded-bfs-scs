# Deterministic Banded BFS for SCS in the Deletion-Only Model

This project implements a deterministic algorithm to compute the shortest common supersequence (SCS) for a set of binary sequences that are known to be derived from a common parent string via deletions. In our model, each input trace is obtained from an unknown parent $$x \in \{0,1\}^n$$ by deleting at most $$t$$ characters. Our algorithm uses a banded breadth-first search (BFS) over a state space defined by pointer tuples (one per trace) to efficiently find the optimal SCS. We also include a known DP merging approach for comparison and a suite of tests (basic, advanced, and timing experiments) to validate correctness, optimality, and performance.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

You will need Python 3 installed. The project also uses the following Python packages:
- `matplotlib`
- `collections` (built-in)
- `itertools` (built-in)

You can install the non-built-in package via pip:
`pip install matplotlib`

### Installing

Clone the repository to your local machine:
`git clone https://github.com/nawafsalman456/banded-bfs-scs.git`

No further installation is required. You can run the Python script directly.

## Running the Tests

To run the automated tests for this system, simply execute the main Python script:
`python scs.py`

This will run:
- **Basic Tests:** These tests verify that for small instances (e.g. with sequences such as `10101` and `1001` for $$t=1$$), the algorithm returns an SCS that is valid (each input is a subsequence) and matches expected output.
- **Advanced Tests:** These tests generate random parent strings and corresponding traces (with exactly or up to $$t$$ deletions) and compare the SCS found by our banded BFS with a known DP merging algorithm. They also check that the deletion invariant is respected.
- **Timing Experiments:** Two timing experiments are conducted:
  - One varying $$n$$ (the length of the parent string) with fixed $$k=3$$ and $$t=5$$.
  - One varying $$k$$ (the number of sequences) with fixed $$n=20$$ and $$t=3$$.
  
  These experiments plot running times to compare our Banded BFS algorithm (which scales nearly linearly in $$n$$ for fixed $$k$$ and small $$t$$) with the DP merging approach (which exhibits quadratic or exponential growth).

### End-to-End Test Example

After running the tests, you should see outputs confirming:
- All basic tests pass (e.g., `["10101", "1001"]` yields SCS `"10101"`).
- Advanced tests report optimality by comparing SCS lengths with the DP merging method.
- Timing experiments display plots saved as `scs_calc_times_vs_n` and `scs_calc_times_vs_k`.

## Deployment

This project is designed as a research prototype. To deploy it on a live system or integrate it into a larger codebase, simply import the provided functions (e.g. `calc_scs_bfs`, `calc_scs_dp`) from the module. You may also adjust the parameters (such as the deletion limit $$t$$ and the banding factor) for your particular application.


## Authors

* **Nawaf Salman**
* **Nahel Awidat**


## Timing Results

Three timing experiments were performed:

1. **Experiment 1: Varying $$n$$ (Length of Parent String) with Fixed $$k=3$$ and $$t=5$$**

   - **Setup:**  
     Parent length $$n$$ varied from 10 to 1000 (in increments of 10). We fixed $$k=3$$ and $$t=5$$. For each $$n$$, a random parent $$x$$ was generated, and 3 traces were produced by deleting exactly $$t$$ symbols from $$x$$.
     
   - **Results:**  
     The Banded BFS algorithm showed nearly linear growth in running time with $$n$$, while the DP merging approach exhibited quadratic growth.
     
   - **Explanation:**  
     The BFS restricts the state space to $$O(n \cdot (t+1)^{k-1})$$ states, which is nearly linear in $$n$$ for fixed $$t$$ and $$k$$. In contrast, the DP merging method performs pairwise merges and explores multiple orderings, leading to quadratic time in $$n$$.

2. **Experiment 2: Varying $$k$$ (Number of Sequences) with Fixed $$n=20$$ and $$t=3$$**

   - **Setup:**  
     With $$n=20$$ and $$t=3$$ fixed, $$k$$ was varied from 3 to 10. A random parent $$x$$ was generated and traces were produced by deleting exactly $$t$$ symbols.
     
   - **Results:**  
     For small $$k$$ (up to about 8), both algorithms were fast. As $$k$$ increased, the DP merging approach’s running time grew dramatically (near factorial), while the Banded BFS algorithm’s running time increased only moderately.
     
   - **Explanation:**  
     The BFS approach exploits the deletion invariant to restrict the state space, keeping the complexity at $$O(n \cdot (t+1)^{k-1})$$. The DP merging approach must try all permutations of sequence merging, resulting in exponential growth in $$k$$.

3. **Experiment 3: Varying $$t$$ (Allowed Deletions) with Fixed $$n=500$$ and $$k=3$$**

   - **Setup:**  
     We fixed the parent length $$n=500$$ and number of sequences $$k=3$$. We varied the deletion parameter $$t$$ from 1 up to 50. For each value of $$t$$, a random parent $$x$$ of length 500 was generated, and three traces were produced by deleting exactly $$t$$     
     symbols from $$x$$.
   
   - **Results:**  
     The Banded BFS is extremely fast for small $$t$$, but as $$t$$ increases, its running time grows significantly. In contrast, the DP merging approach’s running time remains relatively stable around 1–2 seconds for these parameters, occasionally outpacing BFS for 
     very large $$t$$.
   
   - **Explanation:**  
     The Banded BFS approach has a state space of roughly $$O(n \cdot (t+1)^{k-1})$$ for fixed $$k$$, so it scales nearly linearly in $$n$$ but exponentially in $$t$$. Hence, for large $$t$$ (e.g., $$t\approx 50$$), BFS can slow down.  
     The DP merging approach depends more on $$n$$ and $$k$$ and less on $$t$$. Since $$n=500$$ and $$k=3$$ are fixed, its running time stays within a moderate range. Nevertheless, for small $$t$$—the most common scenario in practical applications like DNA analysis—BFS 
     remains significantly faster.
     
---

**Conclusion**

Our deterministic banded BFS algorithm efficiently computes the shortest common supersequence (SCS) in the deletion-only model. Validation via basic and advanced tests confirms both correctness (each trace is a subsequence of the output) and optimality (the output length matches that of the best-known DP method). Timing experiments show that for fixed $$k$$ and small $$t$$, our algorithm scales nearly linearly with the parent string length, and it outperforms the DP merging approach—especially as the number of sequences increases.

## Resources

1. V. I. Levenshtein. “Efficient Reconstruction of Sequences from Their Subsequences or Supersequences.” _J. Combinatorial Theory, Series A_, 93(2):310–332, 2001.
2. G. Navarro. “A Guided Tour to Approximate String Matching.” _ACM Computing Surveys_, 33(1):31–88, 2001.
3. D. Gusfield. _Algorithms on Strings, Trees, and Sequences_. Cambridge University Press, 1997.
4. M. Mitzenmacher. “A Survey of Results for Deletion Channels and Related Synchronization Channels.” _Probability Surveys_, 6:1–33, 2009.
5. J. S. Bae and J. S. Park. “Shortest Common Supersequence Problem: A Survey.” _International Journal of Computers and Their Applications_, 23(4):229–237, 2016.
6. D. Maier. “The Complexity of Some Problems on Subsequences and Supersequences.” _J. ACM_, 25(2):322–336, 1978.
7. R. D. Blum, C. Chalasani, and C. R. N. Rao. “The Complexity of the Shortest Common Supersequence Problem.” _SIAM Journal on Computing_, 24(6):1414–1423, 1995.
8. D. Jiang, L. Wang, and T. Jiang. “A Survey of the Shortest Common Supersequence Problem.” _IEEE/ACM Transactions on Computational Biology and Bioinformatics_, 7(4):758–770, 2010.
