// approx.hpp - approximate matching automata (Hamming and Levenshtein NFAs)
//
// State (i,e) = "i pattern symbols consumed with e errors so far".
//   match         (i,e) --pat[i]-->  (i+1,e)
//   substitution  (i,e) --~pat[i]--> (i+1,e+1)
//   insertion     (i,e) --any-->     (i,e+1)     (extra text symbol, 0<i<m)   [edit only]
//   deletion      (i,e) --eps-->     (i+1,e+1)   (pattern symbol skipped)     [edit only]
#pragma once
#include <climits>
#include "common.hpp"
#include "nfa.hpp"

namespace bio {

enum class ApproxMode { Hamming, Edit };

inline std::vector<Mask> literalMasks(const std::string& p, const Alphabet& al) {
    std::vector<Mask> v;
    for (char c : p)
        if (!std::isspace((unsigned char)c)) v.push_back(al.expand(c));
    return v;
}

inline Nfa buildApproxNfa(const std::vector<Mask>& pat, int k, ApproxMode mode, const Alphabet& al) {
    int m = (int)pat.size();
    if (m == 0) throw ParseError("empty pattern");
    if (k < 0 || k >= m) throw ParseError("need 0 <= k < pattern length");
    Nfa n;
    n.asize = al.size();
    int W = m + 1;
    for (int i = 0; i < W * (k + 1); ++i) n.add();
    auto id = [&](int i, int e) { return e * W + i; };
    n.start = id(0, 0);
    for (int e = 0; e <= k; ++e) {
        for (int i = 0; i <= m; ++i) {
            int q = id(i, e);
            if (i == m) { n.s[q].acc = true; n.s[q].tag = e; continue; }
            n.edge(q, pat[i], id(i + 1, e));
            if (e < k) {
                Mask sub = al.full() & ~pat[i];
                if (sub) n.edge(q, sub, id(i + 1, e + 1));
                if (mode == ApproxMode::Edit) {
                    if (i > 0) n.edge(q, al.full(), id(i, e + 1));
                    n.epsilon(q, id(i + 1, e + 1));  // epsilon edges always go to a larger state id
                }
            }
        }
    }
    return n;
}

struct ApproxHit { size_t start, end; int errors; };

// Streaming simulation of the approximate NFA (an implicit Sigma* prefix: a fresh alignment starts at
// every position). Each active state carries the latest start of any alignment reaching it, so every
// hit gets a span and the minimum number of errors. bestOnly keeps the best non-overlapping hits.
inline std::vector<ApproxHit> approxScan(const Nfa& n, const std::vector<int>& seq, bool bestOnly = true) {
    int N = (int)n.s.size();
    std::vector<long> cur(N, -1), nxt(N, -1);
    auto close = [&](std::vector<long>& v) {
        for (int q = 0; q < N; ++q)
            if (v[q] >= 0) for (int t : n.s[q].eps) v[t] = std::max(v[t], v[q]);
    };
    cur[n.start] = 0;
    close(cur);
    std::vector<ApproxHit> hits;
    for (size_t j = 0; j < seq.size(); ++j) {
        std::fill(nxt.begin(), nxt.end(), -1);
        int a = seq[j];
        if (a >= 0)
            for (int q = 0; q < N; ++q)
                if (cur[q] >= 0)
                    for (auto& e : n.s[q].t)
                        if (e.first >> a & 1u) nxt[e.second] = std::max(nxt[e.second], cur[q]);
        nxt[n.start] = std::max(nxt[n.start], (long)j + 1);
        close(nxt);
        int bestE = INT_MAX; long bestS = -1;
        for (int q = 0; q < N; ++q) {
            if (nxt[q] >= 0 && n.s[q].acc && (size_t)nxt[q] < j + 1) {
                int e = n.s[q].tag;
                if (e < bestE || (e == bestE && nxt[q] > bestS)) { bestE = e; bestS = nxt[q]; }
            }
        }
        if (bestS >= 0) hits.push_back({(size_t)bestS, j + 1, bestE});
        cur.swap(nxt);
    }
    if (!bestOnly) return hits;
    std::vector<ApproxHit> sorted = hits, chosen;
    std::stable_sort(sorted.begin(), sorted.end(), [](const ApproxHit& a, const ApproxHit& b) {
        return a.errors != b.errors ? a.errors < b.errors : a.end < b.end;
    });
    for (auto& h : sorted) {
        bool clash = false;
        for (auto& c : chosen) if (h.start < c.end && c.start < h.end) { clash = true; break; }
        if (!clash) chosen.push_back(h);
    }
    std::sort(chosen.begin(), chosen.end(), [](const ApproxHit& a, const ApproxHit& b) { return a.start < b.start; });
    return chosen;
}

}  // namespace bio
