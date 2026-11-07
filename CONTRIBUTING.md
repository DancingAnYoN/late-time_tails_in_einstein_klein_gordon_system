# Contributing

Thank you for helping develop MMIEF. This is experimental numerical-relativity software.

## Development setup

Use Julia 1.12 (the target version). Run `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'` after modifying the core solver.

The rendering environment is separate. See `README.md` for installation and optional OpenGL smoke checks. The GitHub CI intentionally tests the core without requiring a graphical display.

## Pull requests

- Describe the equations, gauge conventions, and boundary conditions affected.
- Provide a minimal reproducible test and expected versus observed output.
- For scientific changes, include resolution and outer-boundary checks and report constraint norms.
- Do not commit generated `results/`, personal machine paths, or large MP4/CSV outputs.
- Distinguish software regression tests from scientific validation of QNMs or Price tails.

When reporting a bug, include Julia version, operating system, command, stack trace, and a small run configuration.
