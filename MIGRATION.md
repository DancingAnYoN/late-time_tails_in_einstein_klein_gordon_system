# Migration from MMIEF_RK4_GL to MMIEF_Modular

The old `run_simulation.jl` combined numerical evolution and visualization in
one Julia process. The new runner deliberately separates them.

| Previously | Now |
|---|---|
| `include("run_simulation.jl"); main(["--quick", "--no-gl"])` | `include("run_simulation.jl"); main(["--quick"])` |
| `main(["--quick", "--no-movie"])` | First run `main(["--quick"])`, then `julia --project=viz scripts/render_dashboard.jl results/quick` |
| `main(["--quick"])` (dashboard movie inline) | Run solver, then `julia --project=viz scripts/render_dashboard.jl results/quick --movie` |
| `include("pulse_movie.jl"); run_pulse_movie(["--quick"])` | Run solver, then `julia --project=viz scripts/render_movie.jl results/quick` |
| `MMIEFCore.run_evolution(p)` | `MMIEF.run_evolution(p)` |
| `MMIEFCore.save_csv(data,path)` | `MMIEF.save_run(data,run_directory)`; legacy flat `save_csv` still exists |
| `MMIEFCore.fit_qnm(...)` | `MMIEF.fit_qnm(...)` |

Never mix the old `MMIEFCore.jl` with the modular `MMIEF` package in a
production script. The frozen old source lives under `legacy/` for regression
comparisons only. Restart the Julia REPL when switching project versions.

**Note:** `--no-gl`, `--no-movie` and `--show` are not accepted by the new
numerical runner. They are visualization-specific options provided by the
separate scripts.
