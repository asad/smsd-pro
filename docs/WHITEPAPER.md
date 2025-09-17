# SMSD Pro White Paper

## Summary
SMSD Pro implements fast, chemistry‑aware subgraph and MCS search for small molecules. Substructure uses a VF2++‑style state machine; MCS uses a modular product with a bit‑parallel Branch‑and‑Bound Maximum Clique (BBMC) and an optional McGregor grow step. Profiles expose RDKit‑parity behaviours.

### Heritage & Prior Work
- **SMSD (2009)**: Rahman *et al.*, *Journal of Cheminformatics* 1:12. DOI:10.1186/1758-2946-1-12.
- **VF2 family**: Cordella *et al.* (2004). State‑of‑the‑art backtracking with feasibility checks.
- **Maximum Clique**: Tomita *et al.* (2010), San Segundo *et al.* (2013) — colouring bounds and bitsets.
- **McGregor extension**: McGregor (1979) — seed‑and‑grow heuristic improving common subgraphs.

## Chemistry Model
- Atom comparators: element/isotope/charge/valence/aromaticity/ring membership/ring fusion count/stereo.
- Bond comparators: order (strict/loose), stereo (off/defined/exact), ring size tolerance.
- Profiles: `rdkit-strict`, `rdkit-flexible`, `maximal-loose`, `default`.

## Algorithms

### Substructure (VF2++)
1. Pre‑compute compatibility matrix using atom comparators and degree slack.
2. DFS with terminal/frontier counts and feasibility checks on every candidate extension.
3. `connected_only=True` by default; dotted SMARTS fragments allowed when present.
4. Optional induced constraint.

**Pseudo‑code (substructure)**
```
def vf2pp(q, t, C, opt):
    init_state()
    def dfs():
        if mapping_complete(): record; return
        i = pick_next_q(connected=opt.connected_only)
        for j in compatible_targets(i):
            if opt.induced and violates_induced(i, j): continue
            push(i, j)
            dfs()
            pop()
    dfs()
```

### MCS (Modular Product + BBMC)
- Nodes: all (i ∈ mol1, j ∈ mol2) with compatible atoms.
- Edges:  
  **MCIS (induced)** — require `(i,k)` edge iff `(j,l)` edge; connect also on non‑edges.  
  **MCSS/MCCS (non‑induced, connected)** — connect only when both bonds exist *and* are compatible; extra bonds in the target are allowed.
- Maximum clique via BBMC with greedy colour bound and Tomita pivoting.
- Optional **McGregor** grow from the seed clique (skipped for SMARTS queries to avoid pattern valence issues).

**Pseudo‑code (MCS)**
```
def mcs(m1, m2, C, opt):
    P = build_product(m1, m2, C, induced=(opt.type=="MCIS"))
    best = []
    def bbmc(R, Pset, X):
        if upper_bound(Pset) + |R| <= |best|: return
        if not Pset and not X: best = R; return
        u = pivot(Pset ∪ X)
        for v in Pset \ N(u):
            bbmc(R ∪ {v}, Pset ∩ N(v), X ∩ N(v))
            Pset.remove(v); X.add(v)
    bbmc(∅, V(P), ∅)
    mapping = to_mapping(best)
    if opt.type == "MCCS": mapping = largest_connected_in_query(mapping)
    if opt.use_mcgregor_extend and not query_is_smarts: mapping = mcgregor_extend(mapping)
    return mapping
```

## Why it’s fast
- Heavy pruning from chemistry‑aware compatibility and ring size constraints.
- Strong colouring bounds in BBMC.
- Weisfeiler–Lehman (1‑WL) invariants in McGregor provide cheap feasibility filters.
- Caching for recursive SMARTS predicates.

## Future Work
- Parallel BBMC (per‑branch) and improved colour sorting.
- Higher‑order WL (k‑WL) or GNN embeddings as admissible upper bounds.
- Ring system isomorphism constraints (collapsed ring graphs) for speed on PAH‑rich pairs.
- SIMD bitsets and Rust C‑extensions for further acceleration.
