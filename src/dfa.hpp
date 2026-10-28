// dfa.hpp - complete DFA: subset construction, minimisation, boolean closure operations
#pragma once
#include <functional>
#include <map>
#include <queue>
#include <set>
#include "nfa.hpp"

namespace bio {

struct Dfa {
    int n = 0, k = 0, start = 0;
    std::vector<int> d;     // flat transition table, d[q*k + a]
    std::vector<char> acc;

    int next(int q, int a) const { return d[(size_t)q * k + a]; }

    // Subset (powerset) construction. The result is complete: the empty set is the dead state.
    static Dfa fromNfa(const Nfa& nfa) {
        Dfa r;
        r.k = nfa.asize;
        std::map<std::vector<int>, int> id;
        std::vector<std::vector<int>> sets;
        auto get = [&](const std::vector<int>& v) -> int {
            auto it = id.find(v);
            if (it != id.end()) return it->second;
            int i = (int)sets.size();
            id[v] = i;
            sets.push_back(v);
            return i;
        };
        r.start = get(nfa.closure({nfa.start}));
        for (size_t i = 0; i < sets.size(); ++i) {
            std::vector<int> cur = sets[i];
            r.acc.push_back(nfa.accepting(cur));
            for (int a = 0; a < r.k; ++a) r.d.push_back(get(nfa.step(cur, a)));
        }
        r.n = (int)sets.size();
        return r;
    }

    bool accepts(const std::vector<int>& w) const {
        int q = start;
        for (int a : w) {
            if (a < 0) return false;
            q = next(q, a);
        }
        return acc[q];
    }

    // Moore partition refinement.
    Dfa minimized() const {
        std::vector<int> cls(n);
        for (int q = 0; q < n; ++q) cls[q] = acc[q] ? 1 : 0;
        size_t count = std::set<int>(cls.begin(), cls.end()).size();
        while (true) {
            std::map<std::vector<int>, int> sig;
            std::vector<int> nc(n);
            for (int q = 0; q < n; ++q) {
                std::vector<int> v;
                v.push_back(cls[q]);
                for (int a = 0; a < k; ++a) v.push_back(cls[next(q, a)]);
                auto it = sig.find(v);
                if (it == sig.end()) it = sig.emplace(v, (int)sig.size()).first;
                nc[q] = it->second;
            }
            cls = nc;
            if (sig.size() == count) break;
            count = sig.size();
        }
        Dfa r;
        r.k = k; r.n = (int)count;
        r.d.assign((size_t)r.n * k, 0);
        r.acc.assign(r.n, 0);
        for (int q = 0; q < n; ++q) {
            int c = cls[q];
            r.acc[c] = acc[q];
            for (int a = 0; a < k; ++a) r.d[(size_t)c * k + a] = cls[next(q, a)];
        }
        r.start = cls[start];
        return r;
    }

    Dfa complement() const {
        Dfa r = *this;
        for (auto& x : r.acc) x = !x;
        return r;
    }

    // Product construction; f decides acceptance from the two component verdicts.
    static Dfa product(const Dfa& A, const Dfa& B, std::function<bool(bool, bool)> f) {
        Dfa r;
        r.k = A.k;
        std::map<std::pair<int, int>, int> id;
        std::vector<std::pair<int, int>> st;
        auto get = [&](int a, int b) -> int {
            auto it = id.find({a, b});
            if (it != id.end()) return it->second;
            int i = (int)st.size();
            id[{a, b}] = i;
            st.push_back({a, b});
            return i;
        };
        r.start = get(A.start, B.start);
        for (size_t i = 0; i < st.size(); ++i) {
            auto p = st[i];
            r.acc.push_back(f(A.acc[p.first], B.acc[p.second]));
            for (int a = 0; a < r.k; ++a) r.d.push_back(get(A.next(p.first, a), B.next(p.second, a)));
        }
        r.n = (int)st.size();
        return r;
    }
    static Dfa intersect(const Dfa& A, const Dfa& B) { return product(A, B, [](bool x, bool y) { return x && y; }); }
    static Dfa unite(const Dfa& A, const Dfa& B) { return product(A, B, [](bool x, bool y) { return x || y; }); }
    static Dfa minus(const Dfa& A, const Dfa& B) { return product(A, B, [](bool x, bool y) { return x && !y; }); }

    // States from which an accepting state is reachable (everything else is "dead").
    std::vector<char> live() const {
        std::vector<std::vector<int>> rev(n);
        for (int q = 0; q < n; ++q)
            for (int a = 0; a < k; ++a) rev[next(q, a)].push_back(q);
        std::vector<char> lv(n, 0);
        std::queue<int> Q;
        for (int q = 0; q < n; ++q) if (acc[q]) { lv[q] = 1; Q.push(q); }
        while (!Q.empty()) {
            int u = Q.front(); Q.pop();
            for (int p : rev[u]) if (!lv[p]) { lv[p] = 1; Q.push(p); }
        }
        return lv;
    }

    // Shortest accepted word (BFS). Returns false if the language is empty.
    bool shortest(std::vector<int>& w) const {
        std::vector<int> par(n, -1), sym(n, -1);
        std::vector<char> seen(n, 0);
        std::queue<int> Q;
        Q.push(start); seen[start] = 1;
        int found = -1;
        while (!Q.empty()) {
            int u = Q.front(); Q.pop();
            if (acc[u]) { found = u; break; }
            for (int a = 0; a < k; ++a) {
                int v = next(u, a);
                if (!seen[v]) { seen[v] = 1; par[v] = u; sym[v] = a; Q.push(v); }
            }
        }
        if (found < 0) return false;
        w.clear();
        for (int x = found; x != start; x = par[x]) w.push_back(sym[x]);
        std::reverse(w.begin(), w.end());
        return true;
    }
    bool isEmpty() const { std::vector<int> w; return !shortest(w); }
    static bool equivalent(const Dfa& A, const Dfa& B) {
        return product(A, B, [](bool x, bool y) { return x != y; }).isEmpty();
    }

    std::string toDot(const Alphabet& al) const {
        auto lv = live();
        std::ostringstream o;
        o << "digraph DFA {\n  rankdir=LR;\n  node [shape=circle];\n  __s [shape=point];\n  __s -> q" << start << ";\n";
        for (int q = 0; q < n; ++q)
            if (lv[q]) o << "  q" << q << (acc[q] ? " [shape=doublecircle]" : "") << ";\n";
        for (int q = 0; q < n; ++q) {
            if (!lv[q]) continue;
            std::map<int, Mask> by;
            for (int a = 0; a < k; ++a) by[next(q, a)] |= 1u << a;
            for (auto& kv : by)
                if (lv[kv.first]) o << "  q" << q << " -> q" << kv.first << " [label=\"" << al.show(kv.second) << "\"];\n";
        }
        o << "}\n  // transitions into the dead state are omitted\n";
        return o.str();
    }
};

}  // namespace bio
