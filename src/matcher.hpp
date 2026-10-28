// matcher.hpp - running DFAs over sequences
#pragma once
#include "dfa.hpp"

namespace bio {

struct Match { size_t start, end; char strand = '+'; };  // half-open [start, end), 0-based

struct MatchOpts {
    bool overlapping = false;   // report the longest match from every start position
    bool anchorStart = false;   // PROSITE '<'
    bool anchorEnd = false;     // PROSITE '>'
    size_t minLen = 1;
};

// Leftmost-longest matching. For each candidate start the DFA is run until it enters a dead
// state (no accepting state reachable) or the sequence ends; the last accepting position wins.
inline std::vector<Match> findMatches(const Dfa& d, const std::vector<int>& seq, const MatchOpts& o = {}) {
    std::vector<char> lv = d.live();
    std::vector<Match> out;
    size_t n = seq.size();
    size_t i = 0;
    while (i <= n) {
        if (o.anchorStart && i > 0) break;
        int q = d.start;
        long best = -1;
        if (d.acc[q] && (!o.anchorEnd || i == n)) best = (long)i;
        for (size_t j = i; j < n; ++j) {
            int a = seq[j];
            if (a < 0) break;
            q = d.next(q, a);
            if (!lv[q]) break;
            if (d.acc[q] && (!o.anchorEnd || j + 1 == n)) best = (long)j + 1;
        }
        if (best >= 0 && (size_t)best >= i + o.minLen) {
            out.push_back({i, (size_t)best, '+'});
            i = o.overlapping ? i + 1 : (size_t)best;
        } else {
            ++i;
        }
    }
    return out;
}

// One-pass streaming search with the Sigma*.L DFA: reports every position where some match ends.
inline std::vector<size_t> streamEnds(const Dfa& d, const std::vector<int>& seq) {
    std::vector<size_t> ends;
    int q = d.start;
    for (size_t j = 0; j < seq.size(); ++j) {
        int a = seq[j];
        if (a < 0) { q = d.start; continue; }
        q = d.next(q, a);
        if (d.acc[q]) ends.push_back(j + 1);
    }
    return ends;
}

}  // namespace bio
