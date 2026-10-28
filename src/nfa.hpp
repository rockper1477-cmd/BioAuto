// nfa.hpp - epsilon-NFA with symbol-set labelled edges, plus bio-specific operations
#pragma once
#include <sstream>
#include <utility>
#include "common.hpp"

namespace bio {

struct Nfa {
    struct State {
        std::vector<std::pair<Mask, int>> t;  // (symbol set, target)
        std::vector<int> eps;                 // epsilon targets
        bool acc = false;
        int tag = 0;  // free payload (approximate matching stores the error count here)
    };
    std::vector<State> s;
    int start = 0;
    int asize = 4;

    int add() { s.emplace_back(); return (int)s.size() - 1; }
    void edge(int f, Mask m, int t) { s[f].t.push_back({m, t}); }
    void epsilon(int f, int t) { s[f].eps.push_back(t); }

    size_t numTransitions() const {
        size_t n = 0;
        for (auto& st : s) n += st.t.size() + st.eps.size();
        return n;
    }

    // epsilon-closure of a set of states (sorted, unique result)
    std::vector<int> closure(const std::vector<int>& seed) const {
        std::vector<char> seen(s.size(), 0);
        std::vector<int> stack, out;
        for (int q : seed)
            if (!seen[q]) { seen[q] = 1; stack.push_back(q); }
        while (!stack.empty()) {
            int q = stack.back(); stack.pop_back();
            out.push_back(q);
            for (int r : s[q].eps)
                if (!seen[r]) { seen[r] = 1; stack.push_back(r); }
        }
        std::sort(out.begin(), out.end());
        return out;
    }
    // delta-hat on a state set for one symbol, followed by closure
    std::vector<int> step(const std::vector<int>& set, int a) const {
        std::vector<int> nxt;
        for (int q : set)
            for (auto& e : s[q].t)
                if (e.first >> a & 1u) nxt.push_back(e.second);
        return closure(nxt);
    }
    bool accepting(const std::vector<int>& set) const {
        for (int q : set) if (s[q].acc) return true;
        return false;
    }
    // Direct NFA simulation (set-of-states), O(|w| * |Q|^2) worst case.
    bool accepts(const std::vector<int>& w) const {
        auto cur = closure({start});
        for (int a : w) {
            if (a < 0) return false;
            cur = step(cur, a);
            if (cur.empty()) return false;
        }
        return accepting(cur);
    }

    // Reversal: L(R) = { reverse(w) : w in L }.  Flip every edge, new start epsilon-jumps
    // to the old accepting states, old start becomes the only accepting state.
    Nfa reversed() const {
        Nfa r;
        r.asize = asize;
        r.s.resize(s.size() + 1);
        int ns = (int)s.size();
        r.start = ns;
        for (int u = 0; u < (int)s.size(); ++u) {
            for (auto& e : s[u].t) r.s[e.second].t.push_back({e.first, u});
            for (int v : s[u].eps) r.s[v].eps.push_back(u);
            if (s[u].acc) r.s[ns].eps.push_back(u);
        }
        r.s[start].acc = true;
        return r;
    }
    // Reverse complement: L = { revcomp(w) : w in L }. This is how a query is searched on
    // the opposite strand without ever rewriting the sequence.
    Nfa revcomp(const Alphabet& al) const {
        Nfa r = reversed();
        for (auto& st : r.s)
            for (auto& e : st.t) e.first = al.complement(e.first);
        return r;
    }
    // Sigma* . L : streaming search automaton (accepts whenever a match ENDS here).
    Nfa unanchored(Mask full) const {
        Nfa r = *this;
        int ns = r.add();
        r.edge(ns, full, ns);
        r.epsilon(ns, start);
        r.start = ns;
        return r;
    }

    std::string toDot(const Alphabet& al) const {
        std::ostringstream o;
        o << "digraph NFA {\n  rankdir=LR;\n  node [shape=circle];\n  __s [shape=point];\n  __s -> q" << start << ";\n";
        for (size_t i = 0; i < s.size(); ++i)
            o << "  q" << i << (s[i].acc ? " [shape=doublecircle]" : "") << ";\n";
        for (size_t i = 0; i < s.size(); ++i) {
            for (auto& e : s[i].t)
                o << "  q" << i << " -> q" << e.second << " [label=\"" << al.show(e.first) << "\"];\n";
            for (int v : s[i].eps)
                o << "  q" << i << " -> q" << v << " [label=\"eps\", style=dashed];\n";
        }
        o << "}\n";
        return o.str();
    }
};

}  // namespace bio
