# SMSD Pro: A Practical Engine for Substructure and MCS

**Core contributions**

1. **Substructure** uses a **VF2++-inspired** backtracker:
   - Static + dynamic ordering (fewest feasible candidates first).
   - Terminal/frontier look‑ahead, degree and ring feasibility checks.
   - Chemistry gates: element/charge/valence, ring-only matching, ring-size comparators (exact/subset/ignore), aromaticity strict vs flexible, double bond stereo (off/defined/exact), optional chirality.
   - RDKit-style uniqueness by default (unique target atom sets).

2. **MCS (exact induced)**:
   - Build the **modular product** of two molecules using atom/bond compatibility.
   - Solve **Maximum Clique** using **bit-parallel BBMC** with Tomita-like pivoting and a sequential colouring bound.
   - Optionally **post-extend** using a chemistry-aware **McGregor** method, with a WL multiset bound and frontier ordering.

3. **SMARTS** support:
   - We interpret stereo slash/backslash on adjacent single bonds and propagate to central `C=C` as a *defined* stereo requirement when needed.
   - Atom and bond constraints are extracted from SMARTS without calling RDKit matching, and merged with engine predicates.

**Selected references**

- Tomita, Sutani, Higashi, Takahashi, Wakatsuki. *A Simple and Faster Branch-and-Bound Algorithm for Finding a Maximum Clique.* WALCOM 2010.
- Prosser. *Exact Algorithms for Maximum Clique: A Computational Study.* Algorithms, 2012.
- San Segundo et al. *An improved bit parallel exact maximum clique algorithm.* Optimisation Letters, 2013.
- McCreesh & Prosser. *Reducing the Branching in a Branch and Bound Algorithm for the Maximum Clique Problem.* CP 2014.
- Batsyn et al. *Improvements to MCS algorithm for the maximum clique problem.* JCO 2014.
- McCreesh, Prosser, Trimble. *Glasgow Subgraph Solver*, ICGT 2020.
- Hoffmann, McCreesh, Reilly. *Between Subgraph Isomorphism and Maximum Common Subgraph.* AAAI 2017.

The software is independent of RDKit engines; RDKit is used purely as a chemistry kernel (mol objects and I/O).
