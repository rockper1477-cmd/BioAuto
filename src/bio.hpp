// bio.hpp - biological layers on top of the automata: FASTA, reverse complement, ORFs, PROSITE
#pragma once
#include <istream>
#include <map>
#include "dfa.hpp"
#include "matcher.hpp"
#include "regex.hpp"

namespace bio {

struct Record { std::string name, seq; };

inline std::vector<Record> readFasta(std::istream& in) {
    std::vector<Record> r;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') r.push_back({line.substr(1), ""});
        else {
            if (r.empty()) r.push_back({"seq", ""});
            r.back().seq += line;
        }
    }
    return r;
}

inline std::string cleanSeq(const std::string& s, const Alphabet& al) {
    std::string o;
    for (char c : s) {
        if (!std::isalpha((unsigned char)c)) continue;
        c = (char)std::toupper((unsigned char)c);
        if (al.kind == Kind::RNA && c == 'T') c = 'U';
        if (al.kind == Kind::DNA && c == 'U') c = 'T';
        o += c;
    }
    return o;
}

// IUPAC-aware reverse complement of a nucleotide string
inline std::string reverseComplement(const std::string& s) {
    bool rna = s.find_first_of("Uu") != std::string::npos;
    std::string r;
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        char c = (char)std::toupper((unsigned char)*it), x;
        switch (c) {
            case 'A': x = rna ? 'U' : 'T'; break;
            case 'T': case 'U': x = 'A'; break;
            case 'C': x = 'G'; break;  case 'G': x = 'C'; break;
            case 'R': x = 'Y'; break;  case 'Y': x = 'R'; break;
            case 'K': x = 'M'; break;  case 'M': x = 'K'; break;
            case 'B': x = 'V'; break;  case 'V': x = 'B'; break;
            case 'D': x = 'H'; break;  case 'H': x = 'D'; break;
            default: x = c;  // S, W, N self-complementary
        }
        r += x;
    }
    return r;
}

// Standard genetic code, codon index = 16*b1 + 4*b2 + b3 with T=0, C=1, A=2, G=3.
inline char translateCodon(const char* c) {
    static const char* table = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    int v[3];
    for (int i = 0; i < 3; ++i) {
        switch (std::toupper((unsigned char)c[i])) {
            case 'T': case 'U': v[i] = 0; break;
            case 'C': v[i] = 1; break;
            case 'A': v[i] = 2; break;
            case 'G': v[i] = 3; break;
            default: return 'X';
        }
    }
    return table[16 * v[0] + 4 * v[1] + v[2]];
}
inline std::string translate(const std::string& s) {
    std::string p;
    for (size_t i = 0; i + 3 <= s.size(); i += 3) {
        char a = translateCodon(s.c_str() + i);
        if (a == '*') break;
        p += a;
    }
    return p;
}

// The ORF language ATG (non-stop codon)* (TAA|TAG|TGA) is regular. We compile it with our own query compiler;
// "non-stop codon" is spelled out as a union of four codon shapes.
inline const Dfa& orfDfa() {
    static Dfa d = Dfa::fromNfa(compileRegex("ATG([ACG]NN|T[CT]N|TA[CT]|TG[CGT])*(TAA|TAG|TGA)",
                                             Alphabet::make(Kind::DNA))).minimized();
    return d;
}

struct Orf { size_t start, end; char strand; int frame; std::string protein; };

inline std::vector<Orf> findOrfs(const std::string& dna, size_t minCodons) {
    Alphabet al = Alphabet::make(Kind::DNA);
    const Dfa& d = orfDfa();
    size_t n = dna.size();
    std::vector<Orf> res;
    auto scan = [&](const std::string& s, char strand) {
        MatchOpts o; o.overlapping = true;
        size_t lastEnd[3] = {0, 0, 0};
        for (auto& m : findMatches(d, al.encode(s), o)) {
            int fr = (int)(m.start % 3);
            if (m.start < lastEnd[fr]) continue;  // nested in an outer ORF of the same frame
            lastEnd[fr] = m.end;
            if ((m.end - m.start) / 3 < minCodons) continue;
            Orf r;
            r.strand = strand;
            r.frame = fr + 1;
            if (strand == '+') { r.start = m.start; r.end = m.end; }
            else { r.start = n - m.end; r.end = n - m.start; }
            r.protein = translate(s.substr(m.start, m.end - m.start));
            res.push_back(r);
        }
    };
    scan(dna, '+');
    scan(reverseComplement(dna), '-');
    std::sort(res.begin(), res.end(), [](const Orf& a, const Orf& b) { return a.start < b.start; });
    return res;
}

struct PrositePattern { std::string regex; bool nterm = false, cterm = false; };

// PROSITE syntax -> query syntax:  x -> .   {AB} -> [^AB]   e(n,m) -> {n,m}   '<' '>' anchors
inline PrositePattern prositeToRegex(std::string p) {
    PrositePattern r;
    p.erase(std::remove_if(p.begin(), p.end(), [](unsigned char c) { return std::isspace(c); }), p.end());
    if (!p.empty() && p.back() == '.') p.pop_back();
    size_t from = 0;
    while (from <= p.size()) {
        size_t dash = p.find('-', from);
        std::string el = p.substr(from, dash == std::string::npos ? std::string::npos : dash - from);
        from = dash == std::string::npos ? p.size() + 1 : dash + 1;
        if (el.empty()) continue;
        if (el[0] == '<') { r.nterm = true; el.erase(0, 1); }
        if (!el.empty() && el.back() == '>') { r.cterm = true; el.pop_back(); }
        if (el.empty()) continue;
        std::string core;
        size_t i = 0;
        if (el[0] == 'x' || el[0] == 'X') { core = "."; i = 1; }
        else if (el[0] == '[' || el[0] == '{') {
            char cl = el[0] == '[' ? ']' : '}';
            size_t j = el.find(cl);
            if (j == std::string::npos) throw ParseError("unterminated '" + std::string(1, el[0]) + "' in PROSITE pattern");
            std::string body = el.substr(1, j - 1);
            if (body.find_first_of("<>") != std::string::npos) throw ParseError("anchors inside brackets are not supported");
            core = (el[0] == '[' ? "[" : "[^") + body + "]";
            i = j + 1;
        } else if (std::isalpha((unsigned char)el[0])) { core = std::string(1, el[0]); i = 1; }
        else throw ParseError("bad PROSITE element '" + el + "'");
        std::string rest = el.substr(i);
        if (!rest.empty()) {
            if (rest.front() != '(' || rest.back() != ')') throw ParseError("bad repeat in PROSITE element '" + el + "'");
            core += "{" + rest.substr(1, rest.size() - 2) + "}";
        }
        r.regex += core;
    }
    if (r.regex.empty()) throw ParseError("empty PROSITE pattern");
    return r;
}

}  // namespace bio
