# SMSD Pro – Algorithm Overview (teaching notes)

This whitepaper summarises the design choices in SMSD Pro for students of algorithms and cheminformatics.

## Substructure (VF2++)

We implement a VF2++‑style backtracking subgraph isomorphism tailored for molecular graphs:

* **Feasibility checks:** atom type/charge/valence/aromaticity/ring membership; bond order/aromaticity/stereo; ring‑size comparators.
* **Look‑ahead:** terminal set counting (q_term ≤ t_term, q_new ≤ t_new) to prune early.
* **Ordering:** frontier ordering with a degree/ring/aromaticity tie‑breaker.
* **Induced mode:** optional MCIS‑style “no extra edges” constraint during mapping.

For chemistry, we cache ring information and sizes per atom/bond and run a short 1‑WL colour refinement for cheap upper bounds when extending.
Chirality and E/Z are included where requested—by default double‑bond stereo must be **defined** (not necessarily equal), with an **exact** mode available.

## MCS (MCIS/MCCS)

We construct the **modular product** of the two graphs using our chemistry‑aware compatibilities.
Maximum cliques are found with a **bit‑parallel BBMC** variant (Tomita et al.) using a greedy colouring bound and a strong pivot.
The winning clique gives an MCIS mapping; the **MCCS** is the largest connected component of that mapping in the query graph.
Finally, an improved **McGregor** seed‑to‑extend pass greedily extends the mapping using feasibility checks + WL bounds.

## Recursive SMARTS `$()`

We support two forms without altering the engines:

* `$([SMARTS])`: evaluated by running **our own substructure engine** anchored at the candidate atom for the first pattern atom.
* `$name`: named predicates registered in a tiny **VM**; each call is **memoised** by `(mol_id, atom_idx, expr)`.

The atom comparator has a single extra hook `atom_rec_ok(qmol, tmol, qi, tj)` used only when the query is SMARTS containing `$()`.
Everything else (ordering, pruning, McGregor, clique search) remains unchanged and fast.

## Literature (selected)

* Tomita, Sutani, Higashi, Takahashi, Wakatsuki. *A Simple and Faster Branch‑and‑Bound Algorithm for Finding a Maximum Clique.* WALCOM 2010.
* Prosser. *Exact Algorithms for Maximum Clique: A Computational Study.* Algorithms 2012.
* San Segundo et al. *An improved bit parallel exact maximum clique algorithm.* Optimisation Letters 2013.
* McCreesh & Prosser. *Reducing the Branching in a B&B Algorithm for the Maximum Clique Problem.* CP 2014.
* Batsyn et al. *Improvements to MCS algorithm for the maximum clique problem.* JCO 2014.
* Hoffmann, McCreesh, Reilly. *Between Subgraph Isomorphism and Maximum Common Subgraph.* AAAI 2017.
* McCreesh et al. *When Subgraph Isomorphism is Really Hard...* JAIR 2018.

For RDKit usage and comparison we rely on the official documentation and wheels available on PyPI (`rdkit`).

