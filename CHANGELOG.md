# Changelog

## [1.0.0] - 2025-09-16
### Added
- Exact VF2++ substructure engine with chemistry-aware pruning
- MCS via modular product + BBMC maximum clique; optional McGregor extension
- SMARTS-safe WL colours; recursive SMARTS `$()` anchored to current atom
- Visualisation helpers and curated benchmark suite (CSV/JSON/PNGs under `test_output/`)
- Conservative standardisation pipeline (Largest fragment, normalise, reionise, uncharge, canonical tautomer)
- Complete documentation set: README, WHITEPAPER, bench guide, contribution and security policies, CITATION

### Fixed
- Packaging metadata: switched to SPDX license expression (`Apache-2.0`) and removed deprecated license classifiers
- Guarded McGregor on SMARTS queries to avoid RDKit precondition violations

