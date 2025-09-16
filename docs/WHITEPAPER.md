# SMSD Pro: Substructure and MCS via Chemistry‑Aware Search and Maximum Clique

**Version:** 1.0.0  
**Author:** Syed Asad Rahman (BioInception PVT LTD)

## 1. Heritage & Prior Work

SMSD Pro draws inspiration from the original **Small Molecule Subgraph Detector (SMSD)** by Rahman *et al.* (2009), which established a practical toolkit for subgraph detection and MCS in cheminformatics. The present work modernises the architecture for Python, adds a bit‑parallel clique core, recursive SMARTS, and a visual benchmarking suite.

> Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM. *Small Molecule Subgraph Detector (SMSD) toolkit.* **J Cheminformatics** 2009, 1:12. doi:10.1186/1758-2946-1-12.

## 2. Problem Setting

Given molecular graphs \(G=(V,E)\), \(H=(W,F)\) with atom/bond labels, we need:

- **Substructure**: an injective mapping \(f: V 	o W\) preserving labels and adjacencies (induced or not).  
- **MCS (MCIS/MCCS)**: the largest common labelled subgraph between two molecules.

Chemistry imposes label constraints: atom type/charge/isotope/aromaticity/ring context/stereo; bond order/aromaticity/ring size; optional chirality.

## 3. Substructure Engine (VF2++‑style)

We implement a backtracking state machine with frontier ordering and look‑ahead:

1. Precompute compatibility between query atom \(i\) and target atom \(j\) using strict chemistry checks and degree slack.
2. At depth \(d\), select next query atom \(i\) from the frontier (fewest candidates; break ties by degree/ring/aromaticity).
3. Iterate candidate targets \(j\) that satisfy bond‑wise feasibility against already‑mapped neighbours, and advance.
4. Optional induced constraint forbids extra bonds in the target between mapped neighbours.

**Recursive SMARTS** `$()` is anchored to the current atom (first atom in the recursive pattern) and evaluated via a cached predicate VM. This enables patterns like `[C;$(C(=O)O)]` to behave as expected without altering VF2’s core.

### Pseudo‑code (Substructure)

```
function SUBSTRUCTURE_ALL(q, t, chem, opt):
  st ← init_state(q, t, chem)
  results ← ∅ ; seen ← ∅
  def rec():
    if timeout(): return
    if st.finished():
      m ← mapping(st)
      if chem.complete_rings_only and not rings_complete(q, t, m): return
      key ← dedupe_key(m, opt.uniquify_mode)
      if key ∉ seen: results ← results ∪ {m}; seen ← seen ∪ {key}
      return
    i ← pick_next_query_atom(st, opt.connected_only)
    for j in candidate_targets(st, i):
      if opt.induced and not induced_ok(q,t,st,i,j): continue
      st.add(i,j); rec(); st.backtrack()
  rec()
  return results
```

## 4. MCS via Modular Product + Maximum Clique

We construct the **modular product** \(G oxtimes H\): a vertex \((i,j)\) exists when atom labels are compatible; two such vertices are adjacent when the corresponding bonds are both present and compatible **or** both absent (for MCIS). An MCS is a **maximum clique** in this product.

We solve maximum clique with a **bit‑parallel branch‑and‑bound** (BBMC): greedy sequential colouring provides a tight bound; pivoting prunes the candidate set.

Optionally we apply a **McGregor‑style extension**: starting from a maximum clique seed, greedily grow a feasible mapping along the frontier. We **skip** this step for SMARTS queries to avoid undefined valence/H preconditions on RDKit query atoms; the clique result is already valid.

### Pseudo‑code (MCS)

```
function MCS_MAX(m1, m2, chem, opt):
  V, A ← modular_product(m1, m2, chem)
  (best, cliques) ← BBMC_MAX_CLIQUE(A, opt.time_limit_s)
  if best = 0: return ⌀
  M ← mask_to_mapping(V, cliques[0])
  if opt.mcs_type = MCCS: M ← largest_connected_in_query(m1, M)
  if opt.use_mcgregor_extend and not query_is_smarts(m1):
      M ← MCGREGOR_EXTEND(m1, m2, M, chem, opt.extend_time_limit_s)
  return M
```

### Pseudo‑code (BBMC, sketch)

```
function BBMC_MAX_CLIQUE(A):
  best ← 0 ; best_masks ← ∅
  def colour_bound(P):
    colours ← 0
    while P ≠ ∅:
      colours ← colours + 1
      v ← pick_bit(P); X ← ~A[v]
      P ← P \ {v}
      while P ∩ X ≠ ∅:
        u ← pick_bit(P ∩ X)
        X ← X ∧ ~A[u]; P ← P \ {u}
    return colours
  def rec(R, P, X):
    if |R| + colour_bound(P) < best: return
    if P = ∅ and X = ∅:
      update_best(R); return
    u ← pivot(P ∪ X)  # max degree in P
    for v in P \ N(u):
      rec(R ∪ {v}, P ∩ N(v), X ∩ N(v))
      P ← P \ {v}; X ← X ∪ {v}
  rec(∅, V_all, ∅)
  return (best, best_masks)
```

## 5. Chemistry Model

- Atom compatibility: element (wildcards permitted for SMARTS), formal charge, isotope (optional), aromaticity, ring context, optional valence and chirality checks, plus user hook for recursive SMARTS.
- Bond compatibility: order with aromatic equivalence in loose mode, ring co‑membership and optional ring‑size tolerance, double‑bond stereo (defined/exact).

We compute Weisfeiler–Lehman (WL) colours with **SMARTS‑safe** invariants (query atoms return degree instead of valence/H counts) so that reasoning steps never trip RDKit preconditions.

## 6. Correctness and Performance Notes

- BBMC with sequential colouring is a strong exact solver on molecular graph sizes typically encountered in med‑chem.
- Heuristics (frontier ordering, ring/degree pruning, colour bounds) reduce branching substantially.
- Standardisation makes matches more stable across salts/tautomers without being destructive.

## 7. Empirical Bench

`scripts/bench.py` emits per‑case CSV/JSON and substructure/MCS PNGs under `test_output/`. Curated cases include aromatics, stereochemistry, SMARTS with recursive `$()`, and literature “hard” pairs.

## 8. Limitations & Future Work

- Multi‑molecule (disconnected) queries are supported but may need extended induced semantics.
- Chirality and E/Z handling are available but not exhaustively benchmarked in MCS mode.
- Future extensions: ring system isomorphism constraints, parallel BBMC, higher‑order WL invariants, GPU‑accelerated bitsets, and richer recursive predicate libraries.

## 9. References (Selected)

- Rahman SA *et al.* (2009) **SMSD** toolkit. J Cheminformatics 1:12. doi:10.1186/1758-2946-1-12.
- Tomita E, Sutani Y, Higashi T, Takahashi S, Wakatsuki M. *A simple and faster branch-and-bound algorithm for finding a maximum clique.* WALCOM 2010.
- Prosser P. *Exact algorithms for maximum clique: A computational study.* Algorithms 2012.
- San Segundo P, Artíñano E, Rodriguez-Losada D. *An improved bit parallel exact maximum clique algorithm.* Optimization Letters 2013.
- McGregor J. *Backtrack search algorithms and the maximal common subgraph problem.* Software—Practice & Experience 1982.
- McCreesh C, Prosser P. *Reducing the branching in a B&B algorithm for maximum clique.* CP 2014.
- McCreesh C, Delgado-Frias J, Patricio J, Prosser P. *When Subgraph Isomorphism is Really Hard...* JAIR 2018.
