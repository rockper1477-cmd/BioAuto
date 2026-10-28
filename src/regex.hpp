// regex.hpp - query language parser + Thompson construction
//
// Query grammar (alphabet-aware; whitespace ignored):
//   alt    := cat ('|' cat)*
//   cat    := rep*
//   rep    := atom ( '*' | '+' | '?' | '{m}' | '{m,}' | '{m,n}' )*
//   atom   := '(' alt ')' | '[' '^'? letters ']' | '.' | letter
// Letters are IUPAC codes (R,Y,S,W,K,M,B,D,H,V,N) for DNA/RNA, amino acids (plus B,Z,X) for protein.
#pragma once
#include <memory>
#include "nfa.hpp"

namespace bio {

struct Node;
using NodeP = std::shared_ptr<Node>;
struct Node {
    enum Type { Empty, Set, Cat, Alt, Star, Plus, Opt, Rep } type = Empty;
    Mask set = 0;
    NodeP a, b;
    int lo = 0, hi = 0;  // Rep: hi == -1 means unbounded
};

class RegexParser {
public:
    RegexParser(const std::string& s, const Alphabet& al) : src(s), al(al) {}

    NodeP parse() {
        NodeP n = alt();
        skipWs();
        if (pos != src.size()) fail(std::string("unexpected '") + src[pos] + "'");
        return n;
    }

private:
    const std::string src;
    const Alphabet& al;
    size_t pos = 0;

    [[noreturn]] void fail(const std::string& m) const {
        throw ParseError("pattern error at position " + std::to_string(pos) + ": " + m);
    }
    void skipWs() { while (pos < src.size() && std::isspace((unsigned char)src[pos])) ++pos; }
    char peek() { skipWs(); return pos < src.size() ? src[pos] : '\0'; }

    static NodeP mk(Node::Type t, NodeP a = nullptr, NodeP b = nullptr) {
        auto n = std::make_shared<Node>();
        n->type = t; n->a = a; n->b = b;
        return n;
    }
    NodeP mkSet(Mask m) { auto n = mk(Node::Set); n->set = m; return n; }

    NodeP alt() {
        NodeP l = cat();
        while (peek() == '|') { ++pos; NodeP r = cat(); l = mk(Node::Alt, l, r); }
        return l;
    }
    NodeP cat() {
        NodeP l;
        while (true) {
            char c = peek();
            if (c == '\0' || c == '|' || c == ')') break;
            NodeP f = rep();
            l = l ? mk(Node::Cat, l, f) : f;
        }
        return l ? l : mk(Node::Empty);
    }
    int num() {
        skipWs();
        if (pos >= src.size() || !std::isdigit((unsigned char)src[pos])) fail("expected a number");
        int v = 0;
        while (pos < src.size() && std::isdigit((unsigned char)src[pos])) { v = v * 10 + (src[pos] - '0'); ++pos; if (v > 100000) fail("number too large"); }
        return v;
    }
    NodeP rep() {
        NodeP n = atom();
        while (true) {
            char c = peek();
            if (c == '*') { ++pos; n = mk(Node::Star, n); }
            else if (c == '+') { ++pos; n = mk(Node::Plus, n); }
            else if (c == '?') { ++pos; n = mk(Node::Opt, n); }
            else if (c == '{') {
                ++pos;
                int lo = num(), hi = lo;
                if (peek() == ',') { ++pos; hi = (peek() == '}') ? -1 : num(); }
                if (peek() != '}') fail("expected '}'");
                ++pos;
                if (hi != -1 && hi < lo) fail("bad repeat range {m,n} with n < m");
                if (lo > 500 || hi > 500) fail("repeat count too large (max 500)");
                auto r = mk(Node::Rep, n);
                r->lo = lo; r->hi = hi;
                n = r;
            } else break;
        }
        return n;
    }
    NodeP atom() {
        char c = peek();
        if (c == '(') {
            ++pos;
            NodeP n = alt();
            if (peek() != ')') fail("missing ')'");
            ++pos;
            return n;
        }
        if (c == '[') {
            ++pos;
            bool neg = false;
            if (pos < src.size() && src[pos] == '^') { neg = true; ++pos; }
            Mask m = 0;
            while (pos < src.size() && src[pos] != ']') {
                if (!std::isspace((unsigned char)src[pos])) m |= al.expand(src[pos]);
                ++pos;
            }
            if (pos >= src.size()) fail("unterminated '['");
            ++pos;
            if (neg) m = al.full() & ~m;
            if (!m) fail("empty symbol class");
            return mkSet(m);
        }
        if (c == '.') { ++pos; return mkSet(al.full()); }
        if (c == '*' || c == '+' || c == '?' || c == '{' || c == '}' || c == ']')
            fail(std::string("nothing to apply '") + c + "' to");
        ++pos;
        return mkSet(al.expand(c));
    }
};

struct Frag { int s, e; };

inline Frag build(const NodeP& n, Nfa& m) {
    switch (n->type) {
        case Node::Empty: { int s = m.add(), e = m.add(); m.epsilon(s, e); return {s, e}; }
        case Node::Set:   { int s = m.add(), e = m.add(); m.edge(s, n->set, e); return {s, e}; }
        case Node::Cat:   { Frag a = build(n->a, m), b = build(n->b, m); m.epsilon(a.e, b.s); return {a.s, b.e}; }
        case Node::Alt: {
            Frag a = build(n->a, m), b = build(n->b, m);
            int s = m.add(), e = m.add();
            m.epsilon(s, a.s); m.epsilon(s, b.s); m.epsilon(a.e, e); m.epsilon(b.e, e);
            return {s, e};
        }
        case Node::Star: {
            Frag a = build(n->a, m);
            int s = m.add(), e = m.add();
            m.epsilon(s, a.s); m.epsilon(s, e); m.epsilon(a.e, a.s); m.epsilon(a.e, e);
            return {s, e};
        }
        case Node::Plus: {
            Frag a = build(n->a, m);
            int s = m.add(), e = m.add();
            m.epsilon(s, a.s); m.epsilon(a.e, a.s); m.epsilon(a.e, e);
            return {s, e};
        }
        case Node::Opt: {
            Frag a = build(n->a, m);
            int s = m.add(), e = m.add();
            m.epsilon(s, a.s); m.epsilon(s, e); m.epsilon(a.e, e);
            return {s, e};
        }
        case Node::Rep: {
            int s = m.add(), cur = s;
            for (int i = 0; i < n->lo; ++i) { Frag f = build(n->a, m); m.epsilon(cur, f.s); cur = f.e; }
            auto wrap = [&](Node::Type t) { auto w = std::make_shared<Node>(); w->type = t; w->a = n->a; return w; };
            if (n->hi == -1) {
                Frag f = build(wrap(Node::Star), m); m.epsilon(cur, f.s); cur = f.e;
            } else {
                for (int i = 0; i < n->hi - n->lo; ++i) { Frag f = build(wrap(Node::Opt), m); m.epsilon(cur, f.s); cur = f.e; }
            }
            return {s, cur};
        }
    }
    return {0, 0};
}

// Compile a query string into an epsilon-NFA (Thompson construction).
inline Nfa compileRegex(const std::string& pat, const Alphabet& al) {
    RegexParser p(pat, al);
    NodeP ast = p.parse();
    Nfa m;
    m.asize = al.size();
    Frag f = build(ast, m);
    m.start = f.s;
    m.s[f.e].acc = true;
    return m;
}

}  // namespace bio
