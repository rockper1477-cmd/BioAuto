// bioauto - formal-language toolkit for bioinformatics (NFA/DFA compiler, PDA for RNA structure)
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include "approx.hpp"
#include "bio.hpp"
#include "common.hpp"
#include "dfa.hpp"
#include "matcher.hpp"
#include "nfa.hpp"
#include "regex.hpp"
#include "rna.hpp"

using namespace bio;

struct Args {
    std::vector<std::string> pos;
    std::map<std::string, std::string> opt;
    bool has(const std::string& k) const { return opt.count(k) > 0; }
    std::string get(const std::string& k, const std::string& def = "") const { auto it = opt.find(k); return it == opt.end() ? def : it->second; }
    int geti(const std::string& k, int def) const { return has(k) ? std::stoi(opt.at(k)) : def; }
};

static const std::set<std::string> kFlags = {"both", "overlap", "hamming", "no-wobble", "ends"};

static Args parseArgs(const std::vector<std::string>& v, size_t from) {
    Args a;
    for (size_t i = from; i < v.size(); ++i) {
        if (v[i].rfind("--", 0) == 0) {
            std::string k = v[i].substr(2);
            if (kFlags.count(k)) a.opt[k] = "1";
            else {
                if (i + 1 >= v.size()) throw ParseError("option --" + k + " needs a value");
                a.opt[k] = v[++i];
            }
        } else a.pos.push_back(v[i]);
    }
    return a;
}

static void need(const Args& a, size_t n, const char* usage) {
    if (a.pos.size() < n) throw ParseError(std::string("usage: ") + usage);
}

static std::vector<Record> loadRecords(const std::string& arg, const Alphabet& al) {
    std::vector<Record> r;
    std::ifstream f(arg);
    if (f.good()) r = readFasta(f);
    else r.push_back({"query", arg});
    for (auto& x : r) x.seq = cleanSeq(x.seq, al);
    if (r.empty()) throw ParseError("no sequences found in '" + arg + "'");
    return r;
}

static std::string shorten(const std::string& s, size_t max = 40) {
    if (s.size() <= max) return s;
    return s.substr(0, max / 2) + "..." + s.substr(s.size() - (max / 2 - 3));
}

static std::string setStr(const std::vector<int>& v) {
    std::string s = "{";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + std::to_string(v[i]);
    return s + "}";
}

static void printDfaStats(const Nfa& nfa, const Dfa& d, const Dfa& m, const Alphabet& al) {
    std::cout << "  NFA : " << nfa.s.size() << " states, " << nfa.numTransitions() << " transitions\n"
              << "  DFA : " << d.n << " states (subset construction)\n"
              << "  min : " << m.n << " states (Moore minimisation)\n";
    std::vector<int> w;
    if (m.shortest(w)) std::cout << "  shortest member: " << (w.empty() ? "(empty string)" : al.decode(w)) << "\n";
    else std::cout << "  language is empty\n";
}

static void writeOut(const Args& a, const std::string& text) {
    if (a.has("out")) { std::ofstream(a.get("out")) << text; std::cout << "wrote " << a.get("out") << "\n"; }
    else std::cout << text;
}

// ------------------------------------------------------------------ commands
static int cmdCompile(const Args& a) {
    need(a, 1, "compile <pattern> [--alpha dna|rna|protein] [--dot nfa|dfa|min] [--out file]");
    Alphabet al = Alphabet::parse(a.get("alpha", "dna"));
    Nfa nfa = compileRegex(a.pos[0], al);
    Dfa d = Dfa::fromNfa(nfa), m = d.minimized();
    std::cout << "pattern: " << a.pos[0] << "   alphabet: " << al.syms << "\n";
    printDfaStats(nfa, d, m, al);
    if (a.has("dot")) {
        std::string w = a.get("dot");
        writeOut(a, w == "nfa" ? nfa.toDot(al) : w == "dfa" ? d.toDot(al) : m.toDot(al));
    }
    return 0;
}

static int cmdTrace(const Args& a) {
    need(a, 2, "trace <pattern> <string> [--alpha ...]");
    Alphabet al = Alphabet::parse(a.get("alpha", "dna"));
    Nfa nfa = compileRegex(a.pos[0], al);
    Dfa d = Dfa::fromNfa(nfa);
    std::string text = cleanSeq(a.pos[1], al);
    auto w = al.encode(text);
    auto cur = nfa.closure({nfa.start});
    int q = d.start;
    std::cout << "step sym  DFA  NFA active states\n";
    std::cout << "   0   -  q" << q << "   " << setStr(cur) << "\n";
    for (size_t i = 0; i < w.size(); ++i) {
        if (w[i] < 0) { std::cout << "symbol '" << text[i] << "' is outside the alphabet: reject\n"; return 0; }
        cur = nfa.step(cur, w[i]);
        q = d.next(q, w[i]);
        std::cout << "  " << i + 1 << "   " << text[i] << "  q" << q << "   " << setStr(cur) << (cur.empty() ? "  (dead)" : "") << "\n";
    }
    std::cout << (d.acc[q] ? "ACCEPT" : "REJECT") << "  (NFA " << (nfa.accepting(cur) ? "accepts" : "rejects") << ")\n";
    return 0;
}

static int cmdMatch(const Args& a) {
    need(a, 2, "match <pattern> <seq|fasta> [--alpha ..] [--both] [--overlap] [--ends] [--max N]");
    Alphabet al = Alphabet::parse(a.get("alpha", "dna"));
    Nfa nfa = compileRegex(a.pos[0], al);
    Dfa fwd = Dfa::fromNfa(nfa).minimized();
    bool both = a.has("both");
    if (both && !al.nucleic()) throw ParseError("--both only makes sense for DNA/RNA");
    Dfa rev;
    if (both) rev = Dfa::fromNfa(nfa.revcomp(al)).minimized();
    Dfa stream;
    if (a.has("ends")) stream = Dfa::fromNfa(nfa.unanchored(al.full())).minimized();
    MatchOpts o;
    o.overlapping = a.has("overlap");
    size_t maxPrint = (size_t)a.geti("max", 50);
    std::cout << "pattern " << a.pos[0] << " -> minimal DFA with " << fwd.n << " states\n";
    for (auto& rec : loadRecords(a.pos[1], al)) {
        auto enc = al.encode(rec.seq);
        std::vector<Match> ms = findMatches(fwd, enc, o);
        if (both) {
            auto rm = findMatches(rev, enc, o);
            for (auto& m : rm) { m.strand = '-'; ms.push_back(m); }
            std::sort(ms.begin(), ms.end(), [](const Match& x, const Match& y) { return x.start < y.start; });
        }
        std::cout << ">" << rec.name << "  (" << rec.seq.size() << " residues, " << ms.size() << " match(es))\n";
        size_t shown = 0;
        for (auto& m : ms) {
            if (shown++ >= maxPrint) { std::cout << "  ... (use --max to show more)\n"; break; }
            std::string t = rec.seq.substr(m.start, m.end - m.start);
            if (m.strand == '-') t = reverseComplement(t);
            std::cout << "  " << m.strand << "  " << m.start + 1 << ".." << m.end << "  " << shorten(t) << "\n";
        }
        if (a.has("ends")) {
            auto ends = streamEnds(stream, enc);
            std::cout << "  streaming (Sigma*.L DFA, " << stream.n << " states) match end positions:";
            for (size_t e : ends) std::cout << " " << e;
            std::cout << "\n";
        }
    }
    return 0;
}

static int cmdApprox(const Args& a) {
    need(a, 3, "approx <pattern> <k> <seq|fasta> [--hamming] [--alpha ..]");
    Alphabet al = Alphabet::parse(a.get("alpha", "dna"));
    int k = std::stoi(a.pos[1]);
    ApproxMode mode = a.has("hamming") ? ApproxMode::Hamming : ApproxMode::Edit;
    Nfa nfa = buildApproxNfa(literalMasks(a.pos[0], al), k, mode, al);
    Dfa d = Dfa::fromNfa(nfa.unanchored(al.full()));
    std::cout << (mode == ApproxMode::Edit ? "Levenshtein" : "Hamming") << " automaton, k=" << k << ": NFA "
              << nfa.s.size() << " states, streaming DFA " << d.n << " states (" << d.minimized().n << " minimal)\n";
    for (auto& rec : loadRecords(a.pos[2], al)) {
        auto hits = approxScan(nfa, al.encode(rec.seq));
        std::cout << ">" << rec.name << "  (" << hits.size() << " hit(s))\n";
        for (auto& h : hits)
            std::cout << "  " << h.start + 1 << ".." << h.end << "  errors=" << h.errors << "  "
                      << shorten(rec.seq.substr(h.start, h.end - h.start)) << "\n";
    }
    return 0;
}

static int cmdProsite(const Args& a) {
    need(a, 2, "prosite <PROSITE-pattern> <protein|fasta>   e.g. 'N-{P}-[ST]-{P}'");
    Alphabet al = Alphabet::make(Kind::PROTEIN);
    PrositePattern pp = prositeToRegex(a.pos[0]);
    Dfa d = Dfa::fromNfa(compileRegex(pp.regex, al)).minimized();
    std::cout << "PROSITE " << a.pos[0] << "  =>  query '" << pp.regex << "'" << (pp.nterm ? " [N-term]" : "") << (pp.cterm ? " [C-term]" : "")
              << "  (" << d.n << " DFA states)\n";
    MatchOpts o; o.overlapping = true; o.anchorStart = pp.nterm; o.anchorEnd = pp.cterm;
    for (auto& rec : loadRecords(a.pos[1], al)) {
        auto ms = findMatches(d, al.encode(rec.seq), o);
        std::cout << ">" << rec.name << "  (" << ms.size() << " site(s))\n";
        for (auto& m : ms) std::cout << "  " << m.start + 1 << ".." << m.end << "  " << rec.seq.substr(m.start, m.end - m.start) << "\n";
    }
    return 0;
}

static int cmdOrf(const Args& a) {
    need(a, 1, "orf <dna|fasta> [--min codons]");
    Alphabet al = Alphabet::make(Kind::DNA);
    size_t minC = (size_t)a.geti("min", 20);
    std::cout << "ORF automaton: ATG(non-stop codon)*(TAA|TAG|TGA) -> minimal DFA with " << orfDfa().n << " states\n";
    for (auto& rec : loadRecords(a.pos[0], al)) {
        auto orfs = findOrfs(rec.seq, minC);
        std::cout << ">" << rec.name << "  (" << orfs.size() << " ORF(s) >= " << minC << " codons)\n";
        for (auto& o : orfs)
            std::cout << "  " << o.strand << o.frame << "  " << o.start + 1 << ".." << o.end << "  " << (o.end - o.start) / 3 - 1 << " aa  "
                      << shorten(o.protein, 50) << "\n";
    }
    return 0;
}

static int cmdHairpin(const Args& a) {
    need(a, 1, "hairpin <rna|fasta> [--stem-min 4 --stem-max 8 --loop-min 3 --loop-max 8 --no-wobble]");
    Alphabet al = Alphabet::make(Kind::RNA);
    HairpinParams p;
    p.stemMin = a.geti("stem-min", 4); p.stemMax = a.geti("stem-max", 8);
    p.loopMin = a.geti("loop-min", 3); p.loopMax = a.geti("loop-max", 8);
    p.wobble = !a.has("no-wobble");
    Pda pda = buildHairpinPda(p);
    std::cout << "hairpin PDA: " << pda.numStates() << " states, " << pda.numTransitions() << " transitions (stack alphabet {$,A,C,G,U})\n";
    for (auto& rec : loadRecords(a.pos[0], al)) {
        auto hp = findHairpins(rec.seq, p);
        std::cout << ">" << rec.name << "  (" << hp.size() << " stem-loop(s))\n";
        for (auto& h : hp)
            std::cout << "  " << h.start + 1 << ".." << h.end << "  stem=" << h.stem << " loop=" << h.loop << " score=" << h.score << "\n    "
                      << rec.seq.substr(h.start, h.end - h.start) << "\n    " << h.structure << "\n";
    }
    return 0;
}

static int cmdFold(const Args& a) {
    need(a, 1, "fold <rna|fasta>");
    Alphabet al = Alphabet::make(Kind::RNA);
    for (auto& rec : loadRecords(a.pos[0], al)) {
        if (rec.seq.size() > 1500) throw ParseError("sequence too long for the O(n^3) Nussinov DP (max 1500)");
        int pairs;
        std::string db = nussinov(rec.seq, pairs, !a.has("no-wobble"));
        DbResult v = validateDotBracket(rec.seq, db, !a.has("no-wobble"));
        std::cout << ">" << rec.name << "\n" << rec.seq << "\n" << db << "   pairs=" << pairs << "  PDA check: " << (v.ok ? "valid" : v.msg) << "\n";
    }
    return 0;
}

static int cmdDotcheck(const Args& a) {
    need(a, 2, "dotcheck <rna> <dot-bracket>");
    Alphabet al = Alphabet::make(Kind::RNA);
    DbResult r = validateDotBracket(cleanSeq(a.pos[0], al), a.pos[1], !a.has("no-wobble"));
    if (r.ok) std::cout << "valid structure, " << r.pairs << " base pairs" << (r.pseudoknot ? " (pseudoknotted: needs >1 stack, not context-free)" : " (nested: context-free, one stack suffices)") << "\n";
    else std::cout << "INVALID: " << r.msg << "\n";
    return 0;
}

static int cmdCompare(const Args& a) {
    need(a, 2, "compare <pattern1> <pattern2> [--alpha ..]");
    Alphabet al = Alphabet::parse(a.get("alpha", "dna"));
    Dfa A = Dfa::fromNfa(compileRegex(a.pos[0], al)).minimized(), B = Dfa::fromNfa(compileRegex(a.pos[1], al)).minimized();
    std::vector<int> w;
    auto witness = [&](const Dfa& d) { return d.shortest(w) ? (w.empty() ? std::string("(empty)") : al.decode(w)) : std::string("-"); };
    std::cout << "L1: " << A.n << " min-DFA states, L2: " << B.n << " min-DFA states\n";
    std::cout << "L1 == L2 : " << (Dfa::equivalent(A, B) ? "yes" : "no") << "\n";
    Dfa i = Dfa::intersect(A, B), d12 = Dfa::minus(A, B), d21 = Dfa::minus(B, A);
    std::cout << "L1 & L2 empty : " << (i.isEmpty() ? "yes (disjoint)" : "no, e.g. " + witness(i)) << "\n";
    std::cout << "L1 subset of L2 : " << (d12.isEmpty() ? "yes" : "no, " + witness(d12) + " is in L1 only") << "\n";
    std::cout << "L2 subset of L1 : " << (d21.isEmpty() ? "yes" : "no, " + witness(d21) + " is in L2 only") << "\n";
    return 0;
}

static void help() {
    std::cout <<
R"(bioauto - automata & formal languages for sequence analysis
  compile  <pattern> [--alpha dna|rna|protein] [--dot nfa|dfa|min] [--out f]   build NFA/DFA, show sizes
  trace    <pattern> <string>                                                  step-by-step NFA/DFA simulation
  match    <pattern> <seq|fasta> [--both] [--overlap] [--ends] [--max N]       pattern search (IUPAC regex)
  approx   <pattern> <k> <seq|fasta> [--hamming]                               <=k mismatches / edits
  prosite  <PROSITE pattern> <protein|fasta>                                   protein motif search
  orf      <dna|fasta> [--min codons]                                          ORFs on both strands
  hairpin  <rna|fasta> [--stem-min --stem-max --loop-min --loop-max --no-wobble]  PDA stem-loop recognition
  fold     <rna|fasta>                                                         Nussinov structure + PDA validation
  dotcheck <rna> <dot-bracket>                                                 validate a structure string
  compare  <p1> <p2> [--alpha ..]                                              language equivalence / inclusion
  repl                                                                         interactive shell
Query syntax: letters (IUPAC: R Y S W K M B D H V N), [ACG], [^T], . | ( ) * + ? {m,n}
Quote patterns in your shell, e.g.  bioauto match 'GA(A|T)TC' ACGGAATTC
)";
}

static int dispatch(const std::vector<std::string>& t) {
    if (t.empty()) return 0;
    const std::string& c = t[0];
    Args a = parseArgs(t, 1);
    if (c == "compile") return cmdCompile(a);
    if (c == "trace") return cmdTrace(a);
    if (c == "match") return cmdMatch(a);
    if (c == "approx") return cmdApprox(a);
    if (c == "prosite") return cmdProsite(a);
    if (c == "orf") return cmdOrf(a);
    if (c == "hairpin") return cmdHairpin(a);
    if (c == "fold") return cmdFold(a);
    if (c == "dotcheck") return cmdDotcheck(a);
    if (c == "compare") return cmdCompare(a);
    if (c == "help" || c == "--help" || c == "-h") { help(); return 0; }
    throw ParseError("unknown command '" + c + "' (try: help)");
}

int main(int argc, char** argv) {
    std::vector<std::string> t(argv + 1, argv + argc);
    if (t.empty()) { help(); return 0; }
    if (t[0] == "repl") {
        std::cout << "bioauto REPL - type 'help' or 'quit'\n";
        std::string line;
        while (std::cout << "bioauto> " && std::getline(std::cin, line)) {
            std::istringstream is(line);
            std::vector<std::string> tok;
            for (std::string x; is >> x;) tok.push_back(x);
            if (tok.empty()) continue;
            if (tok[0] == "quit" || tok[0] == "exit") break;
            try { dispatch(tok); } catch (const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; }
        }
        return 0;
    }
    try { return dispatch(t); }
    catch (const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 1; }
}
