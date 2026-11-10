# MMIEF agent knowledge pack — installation and use

This package contains **project context** plus **three Agent Skills**.
It is designed for coding agents such as Codex, with self-contained Markdown
references that can also be supplied to other LLMs.

## Install into your existing project

Copy the following **contents**, preserving directory paths, into the root of
`late-time_tails_in_einstein_klein_gordon_system`:

- `AGENTS.md`
- `PROJECT_CONTEXT.md`
- `.agents/skills/mmief-ekg/`
- `.agents/skills/spherical-numerical-relativity/`
- `.agents/skills/black-hole-perturbations/`

Add `SKILLS_README.md` if you want the installation guide in the repository.
No existing Julia source files need to be replaced. Do not copy an archive's
`.git/` directory into an existing clone.

## What to use for each purpose

| Resource | Purpose |
|---|---|
| `PROJECT_CONTEXT.md` | Paste into a new conversation for comprehensive context |
| `AGENTS.md` | Permanent repository instructions to coding agents |
| `mmief-ekg` | Specific code, gauge, equations, outputs, and debugging |
| `spherical-numerical-relativity` | General ADM spherical-reduction and numerical PDE skill |
| `black-hole-perturbations` | Scalar effective potentials, QNMs, tails, and fitting workflows |

To invoke explicitly in Codex, reference the skill by name (for example,
`$mmief-ekg`), or describe a matching task; supported agent versions may
activate the appropriate skill automatically.

## Skill discovery locations

- **Codex (repository):** `.agents/skills/<skill-name>/SKILL.md`.
- **Codex (user-level reusable):** `~/.agents/skills/<skill-name>/SKILL.md`.
- **VS Code Copilot:** project `.github/skills/<skill-name>/SKILL.md`.
- **Claude:** project `.claude/skills/<skill-name>/SKILL.md`.

When copying to another tool, copy the **whole skill folder** including its
`references/` directory; do not copy only `SKILL.md`. Discovery rules differ
by client/version: placing the file in the correct directory enables discovery,
not a guarantee that it will be invoked for every prompt.

## Quick verification

With Python available, from the project root:

```powershell
python scripts/validate_agent_skills.py
```

This bundled check validates filenames, required YAML keys, and internal
Markdown links. It is a lightweight structural check, **not** an official
Agent Skills validator or a Julia runtime/physics test.

For the reference standard, see https://agentskills.io/specification.
For OpenAI Codex skill discovery, see
https://learn.chatgpt.com/docs/build-skills.

## Update strategy

- Keep `AGENTS.md` short and stable.
- Update `PROJECT_CONTEXT.md` after significant equations, numerical approach,
  data format or verification status changes.
- Keep each `SKILL.md` short and procedural; put detailed formulas into
  `references/` and link directly from skill instructions.
- If a skill is meant to work across multiple repositories, copy it to the
  user's global skills directory; eliminate project-specific paths first.
