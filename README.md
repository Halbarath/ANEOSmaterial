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
  monotonic line and reports a non-monotonic one, and *scan*, which finds a
  root on any line and tolerates roundoff at the table edges and on
  plateaus (what the defect-ridden ANEOS tables need);
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
defects producing the documented status codes). The exit status is the
number of failed checks.

## Using it from another repository

Add it once:

    git subtree add --prefix=interpBilinear <url-of-this-repo> main --squash

and later update with

    git subtree pull --prefix=interpBilinear <url-of-this-repo> main --squash

Develop only in this repository, never in a vendored copy, so the copies do
not diverge. Each consumer records the `INTERP_VERSION_*` it requires and
checks it at initialisation, the way EOSlib checks its sublibraries.

## License

GPL-3.0-or-later, see `COPYING`.
