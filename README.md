# interpBilinear

Bilinear interpolation on a rectangular grid, and its inversion along one
axis. This is the interpolation core shared by the tabulated equation of
state libraries ANEOSmaterial, scvheos and reos3,
which are used through EOSlib in pkdgrav3, ballic and tipsytools.

The library is a single header, `interpBilinear.h`, with every function
`static inline`. The consuming libraries vendor it with `git subtree`, so
each one carries its own private copy and nothing clashes at link time when
several of them end up in one binary.

## What it does

* locate a point on the grid (closed intervals, optional roundoff tolerance
  at the edges, see below), evaluate the bilinear interpolant and its two
  partial derivatives in that cell;
* invert the interpolant along one axis at a fixed coordinate on the other
  (T(rho, u), rho(P, T), ...) by bisection on the node index, in two modes
  (see below);
* a nearest-node lookup for integer tables (phase tables);
* second-order central differences at the nodes of a non-uniform axis, for
  building sound-speed tables;
* a monotonicity check of the node lines, for the report at initialisation.

The core is agnostic about what the numbers mean: a library that
interpolates in (log10 rho, log10 T) passes logarithms in and exponentiates
the results, one that interpolates linear values passes them as they are.
It never prints, never exits and never returns a sentinel: everything that
can fail returns a status code (`INTERP_*`, all negative) and leaves its
outputs untouched. The calling library maps the codes to its own failure
contract.

Layout: a table is one contiguous block, x-major, `z[ix*ny + iy]`, with the
inner index `iy` contiguous. In the EOS libraries x is rho and y is T.

## Inversion

Both modes bisect on the node index, O(log n). A bisection whose end nodes
bracket the target keeps a bracketing segment at every step, on any line,
monotonic or not, so it always ends on a root. The modes differ only where
the end nodes do not bracket the target:

* **strict** (`INTERP_INVERT_STRICT`, scvheos and reos3): the target must lie
  in `[f_0, f_n-1]`, otherwise `INTERP_BELOW` or `INTERP_ABOVE`, also where a
  dip or a rise of a non-monotonic line attains it inside;
* **tolerant** (`INTERP_INVERT_TOLERANT`, ANEOSmaterial): a sign bisection
  for lines of either direction; only if the end nodes do not bracket the
  target does it scan for the lowest-index bracketing segment, then apply
  the edge and plateau tolerances, before `INTERP_NOT_ATTAINED`. A target
  equal to `f_0` is node 0.

An infinite target is attained by no segment of a finite line: strict mode
returns `INTERP_ABOVE` for +inf and `INTERP_BELOW` for -inf, tolerant mode
`INTERP_NOT_ATTAINED`. A NaN target is `INTERP_INVALID` in both. Up to 1.0.3
the node tolerance of the tolerant mode accepted an infinite target at node 0
(|inf - f_0| = inf <= nodeTol * inf), so ANEOSmaterial located u = ±inf at
T_min without a fallback.

Neither mode detects a non-monotonic line (that would cost O(n) per call):
there it returns a root, which need not be the physical one. The consumers
report such lines of their tables once at initialisation
(`interpNodeLineMonotonic`). Up to 1.0.3 the header declared
`INTERP_NONMONOTONIC` for strict mode, but the bisection can never end on a
non-bracketing segment, so it was never returned; it is gone.

## Edge tolerance

The table intervals are closed, and a coordinate slightly outside the first
or last node is accepted and clamped onto it if it lies within the band

    tol * |bound| + tolAbs

so that a state on the table edge survives a unit conversion or a
`10^x <-> log10` round trip. The two parts match the two kinds of axis:

* a **linear axis** (ANEOSmaterial: rho, T) carries a relative roundoff, so
  it takes `tol` relative to the bound and `tolAbs = 0`. An absolute part
  would be wrong there: a density axis may start at 1e-25.
* a **log axis** (scvheos, reos3: log10 rho, log10 T) carries the roundoff of
  `log10(10^x)`, which is absolute (~1e-16) as well as relative to |x|, so
  these libraries pass `tolAbs > 0` too (`SCVHEOS_EDGE_TOL_ABS`,
  `REOS3_EDGE_TOL_ABS`). Without it a log bound at or near 0, which the
  arbitrary code-unit shift of log rho can produce, gets no tolerance at all,
  and a table-edge state is rejected for about half of all unit choices.

The functions taking a single `tol` (`interpLocateAxis`, `interpLocate`,
`interpValue`, `interpGradient`, `interpNearest`, `interpInvertY`,
`interpInvertX`) use `tolAbs = 0`. The `*Tol` variants take both parts:

| one `tol` (tolAbs = 0) | both parts |
|---|---|
| `interpLocateAxis` | `interpLocateAxisTol` |
| `interpLocate` | `interpLocateTol` |
| `interpInvertY` | `interpInvertYTol` |
| `interpInvertX` | `interpInvertXTol` |

The single-`tol` functions are the `*Tol` variants with `tolAbs = 0`, bit for
bit (the test suite checks this), so the linear-axis consumers are unaffected.
`tol = tolAbs = 0` disables the tolerance.

The inversions have their own tolerances on the target value (`valTol`,
`edgeTol`, `nodeTol` in `INTERPINVERTOPT`, see the header); the edge
tolerance above applies to the coordinate that is located first.

## Test

    make test

builds and runs `testInterpBilinear`, which needs no table files: it
checks the arithmetic on synthetic grids (a globally bilinear function
reproduced to roundoff, inversion round trips in both modes, constructed
defects producing the documented status codes, the bracket invariant on random
non-monotonic lines, infinite and NaN targets in both modes), the edge tolerance band on
the kinds of axis the consumers use (a linear axis starting at 1e-25, log
axes with a bound at 0, a log axis shifted to 0 by a unit conversion), that
the single-`tol` functions equal the `*Tol` variants with `tolAbs = 0`, and
that the public types still have the layout of `INTERP_ABI_VERSION` (see
below). The exit status is the number of failed checks.

## Using it from another repository

Add it once:

    git subtree add --prefix=interpBilinear <url-of-this-repo> main --squash

and later update with

    git subtree pull --prefix=interpBilinear <url-of-this-repo> main --squash

Develop only in this repository, never in a vendored copy, so the copies do
not diverge.

### Versions and the type layout

Two numbers describe this header:

* `INTERP_VERSION_MAJOR/MINOR/PATCH` (and `INTERP_VERSION_TEXT`), the
  release;
* `INTERP_ABI_VERSION`, the layout of the public types `INTERPGRID`,
  `INTERPCELL` and `INTERPINVERTOPT`. Bump it whenever one of them changes
  (a member added, removed, reordered or retyped), together with the layout
  reference in `testInterpBilinear.c`, which fails until both agree.

The consumers embed `INTERPGRID` in their material structs, and EOSlib
includes the headers of all of them in one translation unit. The include
guard lets only the first vendored copy in there, so a consumer cannot rely
on its own copy being the one in effect. Each consumer header therefore
checks, at compile time and with names carrying its own prefix (the same
unprefixed names in three headers would clash):

    #define ANEOS_INTERP_VERSION_MAJOR_REQUIRED 1
    #define ANEOS_INTERP_VERSION_MINOR_REQUIRED 0
    #define ANEOS_INTERP_VERSION_PATCH_REQUIRED 4
    #define ANEOS_INTERP_ABI_VERSION 1
    #if (INTERP_VERSION_MAJOR * 10000 + INTERP_VERSION_MINOR * 100 + INTERP_VERSION_PATCH) < \
        (ANEOS_INTERP_VERSION_MAJOR_REQUIRED * 10000 + ANEOS_INTERP_VERSION_MINOR_REQUIRED * 100 + ANEOS_INTERP_VERSION_PATCH_REQUIRED)
    #error "ANEOSmaterial: the interpBilinear in effect is older than the required version"
    #endif
    #if !defined(INTERP_ABI_VERSION) || (INTERP_ABI_VERSION != ANEOS_INTERP_ABI_VERSION)
    #error "ANEOSmaterial: the interpBilinear in effect has another type layout"
    #endif

(likewise `SCVHEOS_` and `REOS3_`). The version check is a minimum, the
layout check is exact: a copy with another layout in the same translation
unit stops the compilation instead of making EOSlib and the library read
the same struct differently. After a layout change, pull the new copy into
every consumer and raise its `*_INTERP_ABI_VERSION` (and the required
version) in the same commit.

A consumer that calls a function added in a later release (as scvheos and
reos3 call the `*Tol` variants) must also raise its required version to that
release: in EOSlib the copy in effect may be another consumer's, and an older
copy would otherwise fail with an undeclared function instead of the clear
version error. New functions without a layout change raise the minor version
and leave `INTERP_ABI_VERSION` alone. (The `*Tol` variants came with 1.0.3
without that bump, so a 1.0.3 copy may lack them; the consumers require 1.0.4,
which has them.) Pull the new copy into all three consumers together.

## License

GPL-3.0-or-later, see `COPYING`.