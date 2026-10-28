# bioauto — Formal Automata & Languages for Bioinformatics (pure C++17)

A simulator that **compiles user queries into NFAs/DFAs** for sequence pattern matching and **extends finite
automata** (approximate matching, two-strand search, pushdown automata) for biological sequence analysis,
including **RNA secondary-structure recognition**. No third-party libraries.

## Build & run
```
make            # builds bin/bioauto and bin/tests
make test       # 347 cross-validation checks
./bin/bioauto help
```

## Layout
| file | content |
|---|---|
| `src/common.hpp` | DNA/RNA/protein alphabets, IUPAC codes as bitmask symbol sets, encoding |
| `src/regex.hpp` | query-language parser → AST → **Thompson ε-NFA** |
| `src/nfa.hpp` | NFA, ε-closure, simulation, **reversal, reverse-complement, Σ\*·L**, DOT export |
| `src/dfa.hpp` | **subset construction**, **Moore minimisation**, complement / ∩ / ∪ / difference (**product**), equivalence, shortest witness |
| `src/matcher.hpp` | leftmost-longest search, streaming Σ\*L search |
| `src/approx.hpp` | **Hamming and Levenshtein automata** (k errors), span + error-count reporting |
| `src/rna.hpp` | **pushdown automaton** for stem-loops, Nussinov DP (CYK-style), dot-bracket PDA validator |
| `src/bio.hpp` | FASTA, revcomp, genetic code, **ORF automaton**, PROSITE → query compiler |
| `src/main.cpp` | CLI + REPL |
| `tests/tests.cpp` | brute-force / DP oracles for every construction |

## Theory ↔ biology map
| formal concept | used for |
|---|---|
| regular expression → ε-NFA → DFA → min-DFA | motif queries with IUPAC codes: `TATA[AT]A[AT]`, `GA(A\|T)TC`, `C{2,4}N` |
| NFA reversal + symbol complement | search the opposite strand without rewriting the sequence (`--both`) |
| Σ\*·L automaton | one-pass streaming search over a genome (`--ends`) |
| product construction, complement | `compare`: language equivalence, inclusion, disjointness of two motifs |
| Hamming / Levenshtein NFA | motif search with ≤ k mismatches / indels (`approx`) |
| regular language `ATG(non-stop)*(TAA\|TAG\|TGA)` | ORF finder on 6 frames (`orf`) |
| PROSITE syntax translation | protein motifs with `<` `>` anchors (`prosite`) |
| **pushdown automaton** | hairpins `x·loop·pair(x)ᴿ` — *not regular*; stack stores the stem (`hairpin`) |
| CFG / CYK-style DP (Nussinov) | max-pair nested structure with dot-bracket output (`fold`) |
| multi-stack check | pseudoknots (`[]`) are beyond context-free (`dotcheck`) |

## Examples
```
./bin/bioauto compile 'GA(A|T)TC' --dot min --out m.dot       # dot -Tpng m.dot -o m.png
./bin/bioauto trace   'A(C|G)T' AGT                            # step-by-step NFA/DFA simulation
./bin/bioauto match   'GAATTC' data/dna.fa --both
./bin/bioauto match   'TATA[AT]A[AT]' data/dna.fa --ends
./bin/bioauto approx  TTGACA 1 data/dna.fa [--hamming]
./bin/bioauto prosite 'N-{P}-[ST]-{P}' data/protein.fa
./bin/bioauto orf     data/dna.fa --min 5
./bin/bioauto hairpin data/rna.fa --stem-min 4 --loop-max 8
./bin/bioauto fold    GGGAAAUCCC
./bin/bioauto dotcheck GGAAGGAACCAACC '((..[[..))..]]'
./bin/bioauto compare '[AC]{2}' '(A|C)(A|C)'
./bin/bioauto repl
```
Query syntax: IUPAC letters, `[ACG]`, `[^T]`, `.`, `|`, `( )`, `*`, `+`, `?`, `{m}`, `{m,}`, `{m,n}`.

## Notes / limits
* Sequence characters outside the alphabet (e.g. `N`) are poison symbols: no match can span them.
* Subset construction can blow up exponentially in theory; minimisation and the streaming scan keep typical motifs tiny.
* Per-start matching is O(n·L); `--ends` gives the strict O(n) single pass.
* Hairpin PDA search is exhaustive nondeterminism (cost grows with stem-max); Nussinov is O(n³), capped at 1500 nt.
