// common.hpp - alphabets, IUPAC symbol sets, shared types
#pragma once
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace bio {

// A set of alphabet symbols is a bitmask (bit i  <->  i-th symbol of the alphabet).
// An NFA/DFA edge labelled with a Mask is shorthand for one edge per symbol in the set,
// which keeps IUPAC codes (R = A|G, N = any, ...) cheap.
using Mask = uint32_t;

struct ParseError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

enum class Kind { DNA, RNA, PROTEIN };

struct Alphabet {
    Kind kind = Kind::DNA;
    std::string syms = "ACGT";

    static Alphabet make(Kind k) {
        Alphabet a;
        a.kind = k;
        switch (k) {
            case Kind::DNA: a.syms = "ACGT"; break;
            case Kind::RNA: a.syms = "ACGU"; break;
            case Kind::PROTEIN: a.syms = "ACDEFGHIKLMNPQRSTVWY"; break;
        }
        return a;
    }
    static Alphabet parse(const std::string& n) {
        if (n == "dna") return make(Kind::DNA);
        if (n == "rna") return make(Kind::RNA);
        if (n == "protein" || n == "prot") return make(Kind::PROTEIN);
        throw ParseError("unknown alphabet '" + n + "' (use dna, rna or protein)");
    }

    int size() const { return (int)syms.size(); }
    bool nucleic() const { return kind != Kind::PROTEIN; }
    char symbol(int i) const { return syms[i]; }
    int index(char c) const {
        size_t p = syms.find((char)std::toupper((unsigned char)c));
        return p == std::string::npos ? -1 : (int)p;
    }
    Mask full() const { return (1u << size()) - 1u; }

    // Expand one pattern letter (IUPAC code) into the set of symbols it denotes.
    Mask expand(char c) const {
        char u = (char)std::toupper((unsigned char)c);
        if (kind == Kind::PROTEIN) {
            int i = index(u);
            if (i >= 0) return 1u << i;
            if (u == 'X') return full();
            if (u == 'B') return (1u << index('D')) | (1u << index('N'));
            if (u == 'Z') return (1u << index('E')) | (1u << index('Q'));
            throw ParseError(std::string("unknown amino-acid code '") + c + "'");
        }
        if (kind == Kind::RNA && u == 'U') u = 'T';
        const char* b = nullptr;
        switch (u) {
            case 'A': b = "A"; break;   case 'C': b = "C"; break;
            case 'G': b = "G"; break;   case 'T': b = "T"; break;
            case 'R': b = "AG"; break;  case 'Y': b = "CT"; break;
            case 'S': b = "CG"; break;  case 'W': b = "AT"; break;
            case 'K': b = "GT"; break;  case 'M': b = "AC"; break;
            case 'B': b = "CGT"; break; case 'D': b = "AGT"; break;
            case 'H': b = "ACT"; break; case 'V': b = "ACG"; break;
            case 'N': b = "ACGT"; break;
            default: throw ParseError(std::string("unknown nucleotide code '") + c + "'");
        }
        Mask m = 0;
        for (; *b; ++b) {
            char x = *b;
            if (kind == Kind::RNA && x == 'T') x = 'U';
            m |= 1u << index(x);
        }
        return m;
    }

    // Watson-Crick complement of a symbol set (A<->T/U, C<->G). Index i <-> 3-i.
    Mask complement(Mask m) const {
        Mask r = 0;
        for (int i = 0; i < 4; ++i)
            if (m >> i & 1u) r |= 1u << (3 - i);
        return r;
    }

    std::string show(Mask m) const {
        if (m == 0) return "{}";
        if (m == full()) return ".";
        std::string s;
        for (int i = 0; i < size(); ++i)
            if (m >> i & 1u) s += syms[i];
        return s.size() == 1 ? s : "[" + s + "]";
    }

    // Sequence -> symbol indices. Unknown characters (N, gaps...) become -1 and act as
    // "poison" symbols: no automaton transition exists for them.
    std::vector<int> encode(const std::string& s) const {
        std::vector<int> v;
        v.reserve(s.size());
        for (char c : s) v.push_back(index(c));
        return v;
    }
    std::string decode(const std::vector<int>& w) const {
        std::string s;
        for (int a : w) s += (a >= 0 ? syms[a] : '?');
        return s;
    }
};

}  // namespace bio
