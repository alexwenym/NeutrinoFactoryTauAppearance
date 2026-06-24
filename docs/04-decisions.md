# Decisions log (ADR-lite)

Record physics/design decisions and their rationale as they're made. Each entry:
context → decision → rationale → consequences. Newest first.

## Template
```
## YYYY-MM-DD — <short title>
**Context:** what prompted the decision.
**Decision:** what we chose.
**Rationale:** why.
**Consequences:** what it implies / what it rules out.
```

## 2026-06-18 — Lightweight rewrite instead of nuSQuIDS
**Context:** legacy `../tau_appearance/` study used nuSQuIDS + GENIE on the cluster.
**Decision:** this repo reimplements the chain analytically with a custom C
oscillator (`oschelper`), no heavy dependencies.
**Rationale:** portability, transparency, and tight control over the off-axis flux
and τ-decay modeling that nuSQuIDS doesn't provide out of the box.
**Consequences:** must validate our oscillation/flux against the legacy nuSQuIDS
results; we own the cross-section and decay modeling ourselves.

<!-- Next decision to record: agreed build order (see roadmap proposal). -->
