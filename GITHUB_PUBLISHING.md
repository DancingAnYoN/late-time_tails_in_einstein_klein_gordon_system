# Publishing this repository on GitHub

The repository ZIP contains a complete **source snapshot**. It does not create a remote repository and does not contain a `.git` history, a pinned `Manifest.toml`, or validated numerical outputs.

## Before making it public

1. Review the numerical code and docs and run the Julia tests locally.
2. **Select a software license** (such as MIT, BSD-3-Clause, or GPL-3.0) and add the corresponding `LICENSE` file with the correct copyright owner(s). No license has been assigned automatically.
3. Decide which scientific results, if any, may be shared; generated `results/` are ignored intentionally.
4. Add author and contact details / a `CITATION.cff` only if you wish to publish them.

## Create an empty GitHub repository

Create an **empty** GitHub repository (do not auto-add a README, license, or .gitignore). You may choose a repository name such as `mmief-ekg`.

From PowerShell in the extracted folder, run:

```powershell
git init
git add .
git commit -m "Initial MMIEF source release"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/mmief-ekg.git
git push -u origin main
```

Replace `YOUR-USERNAME` and `mmief-ekg` with the actual GitHub account and repository name. GitHub authentication may use Git Credential Manager, SSH, or `gh auth login`.

If the repository already exists and contains commits, do not force-push; reconcile its history before pushing.

## GitHub Actions

`.github/workflows/ci.yml` runs core tests, a synthetic QNM interpreter test, and a small vacuum evolution on Ubuntu and Windows. The GLMakie visualization tests are intentionally manual because they need an OpenGL-capable graphical environment. **The CI workflow has been configured but not run in this build environment.**

## Keeping a reproducible research record

For any scientific result used in a paper or presentation, archive the exact source commit, complete `config.toml`/`metadata.toml`, Julia version, resolved dependency manifests, run outputs, and convergence diagnostics. The main `Project.toml` is not a substitute for a frozen computational environment.
