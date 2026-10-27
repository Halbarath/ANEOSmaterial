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
  at the edges), evaluate the bilinear interpolant and its two partial
  derivatives in that cell;
* invert the interpolant along one axis at a fixed coordinate on the other
  (T(rho, u), rho(P, T), ...), in two modes: *strict*, which assumes a
  monotonic line and reports a non-monotonic one, and *tolerant*, which
  finds a root on any line (bisection first, a scan only as a fallback) and
  tolerates roundoff at the table edges and on plateaus (what the
  defect-ridden ANEOS tables need);
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

## Test

    make test

builds and runs `testInterpBilinear`, which needs no table files: it
checks the arithmetic on synthetic grids (a globally bilinear function
reproduced to roundoff, inversion round trips in both modes, constructed
defects producing the documented status codes), and that the public types
still have the layout of `INTERP_ABI_VERSION` (see below). The exit status
is the number of failed checks.

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
    #define ANEOS_INTERP_VERSION_PATCH_REQUIRED 3
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

## License

GPL-3.0-or-later, see `COPYING`.