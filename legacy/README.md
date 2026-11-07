# Frozen legacy solver

This is an unmodified copy of the previous monolithic `MMIEFCore.jl` from
`MMIEF_RK4_GL` (RK4 multiplication and diagnostic fixes included). It is kept
only as a regression baseline; the active package lives in `src/`.

Do **not** include this file in ordinary simulations. It defines a separate
`MMIEFCore` module. The `test/regression/test_legacy_equivalence.jl` harness
loads it in a test process and compares outputs against the modular solver.
