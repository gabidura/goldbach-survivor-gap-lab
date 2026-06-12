# Method notes

For each even integer `N`, the program sets `z = floor(sqrt(N))` and searches the complete-sieve interior interval

```text
z < x < N - z
```

for a prime pair

```text
x + (N - x) = N.
```

A successful interior pair is a complete-sieve survivor representation. The interior condition is important: if both coordinates exceed `sqrt(N)` and survive all divisibility tests by primes up to `sqrt(N)`, then both coordinates are prime.

The program supports several modes:

- `first-interior`: exhaustive scan over even `N`, recording first interior prime pairs and record cases.
- `gap-sample`: sampled complete-window gap diagnostics.
- `exact-gap`: exact complete-window gap diagnostics over smaller ranges.
- `single`: single-`N` diagnostic, usually with Miller--Rabin for very large `N`.

The dense backend builds an odd-only prime bitset up to the largest tested `N`. Memory is approximately `stopN / 16` bytes.

The Miller backend uses deterministic 64-bit Miller--Rabin and does not store a dense prime table.
