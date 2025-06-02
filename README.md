# lichao

Implementation of a Li-Chao tree ported from the authors competitive programming library.

Theoretically, a Li-Chao tree should support any function which has the transcending property, but this implementation supports only lines, which are the most common use case.

Currently doesn't support adding custom-width line segments, will be added soon.

Since the performance of Li-Chao trees depends on the size of the domain, it may be preferable to use the Convex hull trick instead.

# Li-Chao trees
Li-Chao trees solve the following problem class in O(log n) time:
```
Given a set S of functions s.t. each pair of functions can intersect at most one point (e.g. lines), efficiently implement
1. Adding a function to S
2. Answering the minimum/maximum value at f(t) for all functions f in S (i.e. \min{f \in S} f(t))
```

Formally, we require the functions to have the *transcending property*, that is, given two functions $f, g \in S$, if $f(t) = g(t)$, then either $f(x) < g(x) \forall x < t$ and $f(x) > g(x) \forall x > t$, or $f(x) > g(x) \forall x < t$ and $f(x) < g(x) \forall x > t$
