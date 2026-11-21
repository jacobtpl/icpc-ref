# icpc-ref

Our team reference for online contests (UCup), forked from
[KACTL](https://github.com/kth-competitive-programming/kactl).
Code is copy-pasted from [content/](./content/) during contests, so there is no PDF and no page limit.

## Principles

- **Copy-paste.** Each file is pasted after [the template](./content/contest/template.cpp).
  If a file `#include`s another file from the reference, paste that one first.
- **Plain-text comments.** A header comment says what the code does, how to call it and
  how to adapt it, and nothing more. No LaTeX.
- **Adaptability first.** Code should be easy to adapt to a different problem, and hard
  to break while adapting it. This matters more than being short.
- **Portability.** Everything compiles with g++ on Apple Silicon (local testing) and on
  Linux x86-64 (the judge). Linux performance is the one that matters.
  x86-only files say so in their header.

## File headers

```
/**
 * Author: ...
 * Source: ...
 * License: ...
 * Description: What it does, in plain text.
 * Usage: how to call it (only when not obvious)
 * Time: O(N log N)
 * Status: how it was tested
 */
```

Write math the way you would in a code comment: `O(N log N)`, `a^b mod c`, `x <= 1e9`, `sum_{d|n} phi(d) = n`.

The `chapter.tex` files are left over from the KACTL PDF; they still hold the math notes
(formulas, bit tricks) but are not built any more.

## Testing

- `make test-compiles` compiles every header in `content/` on its own, after the template
  (`doc/scripts/skip_headers` lists the few that are snippets rather than complete code).
- `make test` builds and runs the stress tests in `stress-tests/`, which compare against
  naive implementations on random cases.
- Set `CXX` to choose the compiler, e.g. `CXX=g++-16 make test-compiles` with Homebrew g++ on macOS.

CI runs both on Linux and on an Apple Silicon macOS runner. The macOS stress-test step is
informational for now: some geometry tests still assume an 80-bit long double.

## License

As usual for competitive programming, the licensing situation is a bit unclear.
Many source files are marked with license (we try to go with
[CC0](https://creativecommons.org/share-your-work/public-domain/cc0/)), but many also aren't.
Presumably good will is to be assumed from other authors, though, and in many cases permission should not be needed since the code is not distributed.
To help trace things back, sources and authors are noted in source files.

Everything in `stress-tests` is implicitly CC0, except reference implementations taken from around the Internet.
