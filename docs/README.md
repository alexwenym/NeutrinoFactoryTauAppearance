# Documentation — NeutrinoFactoryTauAppearance

Working docs for the project. Maintained by Claude as work proceeds; meant to be
skimmable so progress is visible at a glance.

## The one-line question
Can a neutrino factory probe **δ_CP** via **ν_e → ν_τ appearance**, detected by
precisely measuring **off-axis muons from the leptonic decay τ → μ ν_μ ν_τ** in a
far detector?

## Map of these docs

| File | What's in it |
|------|--------------|
| [00-physics.md](00-physics.md) | The physics goal, the measurement strategy, and the signal/background picture |
| [01-code-map.md](01-code-map.md) | Module-by-module map of the current codebase (what each file does today) |
| [02-roadmap.md](02-roadmap.md) | What's missing, planned work, and the proposed build order |
| [03-progress.md](03-progress.md) | Chronological log of what's been done (newest first) |
| [04-decisions.md](04-decisions.md) | Design/physics decisions and their rationale (ADR-lite) |
| [05-open-questions.md](05-open-questions.md) | Unresolved physics & feasibility questions |
| [06-simulation-design.md](06-simulation-design.md) | The Monte Carlo chain, v1 scope, and the MCS feasibility caveat |
| [07-cluster-handoff.md](07-cluster-handoff.md) | **Start here on the cluster** — current state, env setup, and the alouette swap task |

## Conventions
- **Status tags:** ✅ done · 🚧 in progress · ⬜ not started · ❓ open question · 🐛 known bug
- Dates are absolute (YYYY-MM-DD).
- When code changes, update [03-progress.md](03-progress.md) and the relevant section
  of [01-code-map.md](01-code-map.md) / [02-roadmap.md](02-roadmap.md).
