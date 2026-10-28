// Self-checking tests: every automaton construction is cross-validated against brute force / DP oracles.
#include <iostream>
#include <random>
#include "approx.hpp"
#include "bio.hpp"
#include "dfa.hpp"
#include "matcher.hpp"
#include "regex.hpp"
#include "rna.hpp"

using namespace bio;

static int g_pass = 0, g_fail = 0;
#define CHECK(c) do { if (c) ++g_pass; else { ++g_fail; std::cerr << "FAIL line " << __LINE__ << ": " #c "\n"; } } while (0)

static std::vector<std::vector<int>> allWords(int k, int maxLen) {
    std::vector<std::vector<int>> r{{}};
    size_t b = 0;
    for (int len = 1; len <= maxLen; ++len) {
        size_t e = r.size();
        for (size_t i = b; i < e; ++i)
            for (int a = 0; a < k; ++a) { auto w = r[i]; w.push_back(a); r.push_back(w); }
        b = e;
    }
    return r;
}
static std::string randomSeq(std::mt19937& g, const std::string& syms, size_t n) {
    std::string s;
    for (size_t i = 0; i < n; ++i) s += syms[g() % syms.size()];
    return s;
}

int main() {
    Alphabet dna = Alphabet::make(Kind::DNA);
    std::mt19937 rng(12345);

    // 1. parser / semantics
    {
        Nfa n = compileRegex("A(C|G)T", dna);
        CHECK(n.accepts(dna.encode("ACT"))); CHECK(n.accepts(dna.encode("AGT")));
        CHECK(!n.accepts(dna.encode("ATT"))); CHECK(!n.accepts(dna.encode("ACTT")));
        Nfa r = compileRegex("RYN", dna);  // IUPAC
        CHECK(r.accepts(dna.encode("ACG"))); CHECK(r.accepts(dna.encode("GTT"))); CHECK(!r.accepts(dna.encode("CCA")));
        Nfa q = compileRegex("AC{2,3}G", dna);
        CHECK(!q.accepts(dna.encode("ACG"))); CHECK(q.accepts(dna.encode("ACCG"))); CHECK(q.accepts(dna.encode("ACCCG"))); CHECK(!q.accepts(dna.encode("ACCCCG")));
        Nfa neg = compileRegex("[^T]+", dna);
        CHECK(neg.accepts(dna.encode("ACGGA"))); CHECK(!neg.accepts(dna.encode("ACT")));
        bool threw = false;
        try { compileRegex("A(C", dna); } catch (const ParseError&) { threw = true; }
        CHECK(threw);
        threw = false;
        try { compileRegex("AUG", dna); } catch (const ParseError&) { threw = true; }
        CHECK(threw);
        CHECK(compileRegex("AUG", Alphabet::make(Kind::RNA)).accepts(Alphabet::make(Kind::RNA).encode("AUG")));
    }

    // 2. NFA == DFA == minimal DFA on every word up to length 7
    {
        std::vector<std::string> pats = {"A(C|G)*T", "(AC|GT)+N?", "[AG]{2,3}C", "(A|C)*G(A|C){2}", "TATA[AT]A[AT]", "(ACG)*|T+", "N{3}"};
        auto words = allWords(4, 7);
        for (auto& p : pats) {
            Nfa n = compileRegex(p, dna);
            Dfa d = Dfa::fromNfa(n), m = d.minimized();
            CHECK(Dfa::equivalent(d, m)); CHECK(m.n <= d.n);
            bool agree = true;
            for (auto& w : words) { bool x = n.accepts(w); if (x != d.accepts(w) || x != m.accepts(w)) agree = false; }
            CHECK(agree);
        }
    }

    // 3. closure operations vs brute force
    {
        Nfa a = compileRegex("A(C|G)*", dna), b = compileRegex("(A|C)+G?", dna);
        Dfa A = Dfa::fromNfa(a).minimized(), B = Dfa::fromNfa(b).minimized();
        Dfa I = Dfa::intersect(A, B), U = Dfa::unite(A, B), M = Dfa::minus(A, B), C = A.complement();
        bool ok = true;
        for (auto& w : allWords(4, 6)) {
            bool x = A.accepts(w), y = B.accepts(w);
            if (I.accepts(w) != (x && y) || U.accepts(w) != (x || y) || M.accepts(w) != (x && !y) || C.accepts(w) != !x) ok = false;
        }
        CHECK(ok);
        CHECK(Dfa::equivalent(Dfa::fromNfa(compileRegex("(A|C)(A|C)", dna)).minimized(), Dfa::fromNfa(compileRegex("[AC]{2}", dna)).minimized()));
        CHECK(!Dfa::equivalent(A, B));
        CHECK(Dfa::intersect(Dfa::fromNfa(compileRegex("A+", dna)), Dfa::fromNfa(compileRegex("C+", dna))).isEmpty());
    }

    // 4. reverse complement automaton: w in L  <=>  revcomp(w) in L(rc)
    {
        for (const char* p : {"GAATTC", "A(C|G)+T", "[AG]{2}TN", "ACG|TTT+"}) {
            Nfa n = compileRegex(p, dna), rc = n.revcomp(dna);
            bool ok = true;
            for (auto& w : allWords(4, 6)) {
                std::string s = dna.decode(w);
                if (n.accepts(w) != rc.accepts(dna.encode(reverseComplement(s)))) ok = false;
            }
            CHECK(ok);
        }
        // palindromic restriction site EcoRI is its own reverse complement
        Dfa e = Dfa::fromNfa(compileRegex("GAATTC", dna)), er = Dfa::fromNfa(compileRegex("GAATTC", dna).revcomp(dna));
        CHECK(Dfa::equivalent(e, er));
    }

    // 5. matcher and streaming automaton
    {
        Nfa n = compileRegex("GA(A|T)TC", dna);
        Dfa d = Dfa::fromNfa(n).minimized();
        auto enc = dna.encode("CCGAATCAGATTCGGAATC");
        auto ms = findMatches(d, enc);
        CHECK(ms.size() == 3); CHECK(ms[0].start == 2 && ms[0].end == 7);
        std::string patterns[] = {"A(C|G)+T", "TA?T", "(GA)+"};
        for (auto& p : patterns) {
            Nfa nn = compileRegex(p, dna);
            Dfa sd = Dfa::fromNfa(nn.unanchored(dna.full())).minimized();
            for (int rep = 0; rep < 20; ++rep) {
                std::string s = randomSeq(rng, "ACGT", 30);
                auto w = dna.encode(s);
                auto ends = streamEnds(sd, w);
                std::set<size_t> got(ends.begin(), ends.end()), want;
                for (size_t e = 1; e <= w.size(); ++e)
                    for (size_t st = 0; st < e; ++st)
                        if (nn.accepts(std::vector<int>(w.begin() + st, w.begin() + e))) { want.insert(e); break; }
                CHECK(got == want);
            }
        }
        // N in the sequence blocks matches
        CHECK(findMatches(d, dna.encode("GANTC")).empty());
    }

    // 6. approximate matching vs Sellers-style DP oracle
    {
        for (int rep = 0; rep < 60; ++rep) {
            std::string pat = randomSeq(rng, "ACGT", 4 + rng() % 4), text = randomSeq(rng, "ACGT", 40);
            int m = (int)pat.size(), k = 1 + (int)(rng() % 2);
            auto pm = literalMasks(pat, dna);
            auto w = dna.encode(text);
            // edit distance
            Nfa ne = buildApproxNfa(pm, k, ApproxMode::Edit, dna);
            auto hits = approxScan(ne, w, false);
            std::map<size_t, int> got;
            for (auto& h : hits) got[h.end] = h.errors;
            std::vector<std::vector<int>> D(m + 1, std::vector<int>(text.size() + 1, 0));
            for (int i = 1; i <= m; ++i) D[i][0] = i;
            for (size_t j = 1; j <= text.size(); ++j)
                for (int i = 1; i <= m; ++i) {
                    int v = std::min(D[i - 1][j - 1] + (pat[i - 1] == text[j - 1] ? 0 : 1), D[i - 1][j] + 1);
                    if (i < m) v = std::min(v, D[i][j - 1] + 1);
                    D[i][j] = v;
                }
            bool ok = true;
            for (size_t j = 1; j <= text.size(); ++j) {
                bool want = D[m][j] <= k;
                if (want != (got.count(j) > 0)) ok = false;
                else if (want && got[j] != D[m][j]) ok = false;
            }
            CHECK(ok);
            // Hamming
            Nfa nh = buildApproxNfa(pm, k, ApproxMode::Hamming, dna);
            auto hh = approxScan(nh, w, false);
            std::map<size_t, int> gh;
            for (auto& h : hh) gh[h.end] = h.errors;
            ok = true;
            for (size_t j = 1; j <= text.size(); ++j) {
                int mm = -1;
                if ((int)j >= m) { mm = 0; for (int i = 0; i < m; ++i) mm += pat[i] != text[j - m + i]; }
                bool want = mm >= 0 && mm <= k;
                if (want != (gh.count(j) > 0)) ok = false;
                else if (want && gh[j] != mm) ok = false;
            }
            CHECK(ok);
        }
        // k=0 equals exact
        Nfa ex = buildApproxNfa(literalMasks("ACGT", dna), 0, ApproxMode::Edit, dna);
        CHECK(approxScan(ex, dna.encode("TTACGTTT")).size() == 1);
    }

    // 7. hairpin PDA vs brute force
    {
        Alphabet rna = Alphabet::make(Kind::RNA);
        for (bool wob : {true, false}) {
            HairpinParams p; p.stemMin = 2; p.stemMax = 4; p.loopMin = 3; p.loopMax = 5; p.wobble = wob;
            Pda pda = buildHairpinPda(p);
            for (int rep = 0; rep < 25; ++rep) {
                std::string s = randomSeq(rng, "ACGU", 30);
                std::set<std::pair<size_t, size_t>> want, got;
                for (size_t i = 0; i < s.size(); ++i) {
                    for (size_t len : pda.acceptingPrefixes(s, i)) got.insert({i, len});
                    for (int n = p.stemMin; n <= p.stemMax; ++n)
                        for (int l = p.loopMin; l <= p.loopMax; ++l) {
                            size_t len = 2 * n + l;
                            if (i + len > s.size()) continue;
                            bool ok = true;
                            for (int t = 0; t < n; ++t) ok = ok && canPair(s[i + t], s[i + len - 1 - t], wob);
                            if (ok) want.insert({i, len});
                        }
                }
                CHECK(got == want);
            }
        }
        HairpinParams p;
        auto hp = findHairpins("AAGGGGAAAACCCCAA", p);
        CHECK(hp.size() == 1); CHECK(hp[0].stem == 4 && hp[0].loop == 4);
        CHECK(findHairpins("ACACACACACACAC", p).empty());
    }

    // 8. Nussinov and dot-bracket
    {
        int pairs = 0;
        std::string db = nussinov("GGGAAAUCCC", pairs);
        CHECK(pairs == 3); CHECK(db == "(((....)))" || validateDotBracket("GGGAAAUCCC", db).ok);
        for (int rep = 0; rep < 20; ++rep) {
            std::string s = randomSeq(rng, "ACGU", 40);
            std::string d = nussinov(s, pairs);
            DbResult r = validateDotBracket(s, d);
            CHECK(r.ok); CHECK(r.pairs == pairs);
        }
        CHECK(validateDotBracket("GGGAAAUCCC", "(((....)))").ok);
        CHECK(!validateDotBracket("GGGAAAUCCC", "(((....))").ok);
        CHECK(!validateDotBracket("GGGAAAUCCC", "(((....)))x").ok);
        CHECK(!validateDotBracket("AAAAAAAAAA", "(((....)))").ok);   // A-A not a pair
        CHECK(!validateDotBracket("GGGAAUCCC", "(((..)))").ok);       // length mismatch
        DbResult pk = validateDotBracket("GGAAGGAACCAACC", "((..[[..))..]]", true, 2);
        CHECK(pk.ok); CHECK(pk.pseudoknot);
    }

    // 9. ORFs, translation, PROSITE
    {
        auto orfs = findOrfs("CCATGAAATAGCC", 1);
        CHECK(orfs.size() == 1); CHECK(orfs[0].start == 2 && orfs[0].end == 11); CHECK(orfs[0].protein == "MK");
        // reverse strand: revcomp of ATGAAATAG is CTATTTCAT
        auto rv = findOrfs("CCCTATTTCATCC", 1);
        CHECK(rv.size() == 1); CHECK(rv[0].strand == '-'); CHECK(rv[0].protein == "MK");
        CHECK(findOrfs("ATGAAACCC", 1).empty());               // no stop codon
        CHECK(translate("ATGGCCTGGTAA") == "MAW");
        Alphabet pr = Alphabet::make(Kind::PROTEIN);
        PrositePattern pp = prositeToRegex("N-{P}-[ST]-{P}.");
        CHECK(pp.regex == "N[^P][ST][^P]");
        Dfa d = Dfa::fromNfa(compileRegex(pp.regex, pr)).minimized();
        MatchOpts o; o.overlapping = true;
        auto ms = findMatches(d, pr.encode("MNASAKNPSA"), o);
        CHECK(ms.size() == 1); CHECK(ms[0].start == 1 && ms[0].end == 5);
        PrositePattern an = prositeToRegex("<M-x(2,3)-K>");
        CHECK(an.nterm && an.cterm && an.regex == "M.{2,3}K");
        Dfa ad = Dfa::fromNfa(compileRegex(an.regex, pr)).minimized();
        MatchOpts ao; ao.overlapping = true; ao.anchorStart = true; ao.anchorEnd = true;
        CHECK(findMatches(ad, pr.encode("MAAK"), ao).size() == 1);
        CHECK(findMatches(ad, pr.encode("AMAAK"), ao).empty());
        CHECK(findMatches(ad, pr.encode("MAAKA"), ao).empty());
    }

    std::cout << g_pass << " checks passed, " << g_fail << " failed\n";
    return g_fail ? 1 : 0;
}
