# Shared research instructions (Claude Code and Codex)

## Always-on principles

- **Current task first.** The user's current request and the active conversation define what to solve, the scope, and when to stop. Project documents supply standards and evidence; they are not additional assignments or instructions to reopen historical work.
- **Scientific goal.** Study reliable, high-precision robotic contact-force control in dynamic, spatially nonuniform environments (including chest/abdominal ultrasound scanning). The prior B0 proposal is a revisable research starting point, not a mandatory solution or a reason to change the user's current subproblem.
- **Research standard.** Work comes first; a rigorous, coherent Method and a supported Story also matter. Start from the physical problem, uncover useful mathematical structure and hidden conditions, use or adapt mature theory when appropriate, and design new mechanisms only when justified. Do not choose algorithms for novelty or complexity.
- **Evidence standard.** Separate established facts, hypotheses, candidates, implementations and tested cases. Link theoretical claims to assumptions and derivations, and empirical claims to reproducible evidence and actual control benefit. Investigate disagreements rather than tuning simulations to fit expectations. Do not extrapolate local results to an entire research direction, or use locality to protect a failing idea indefinitely.
- **Efficient work.** Use the least work that can answer the current scientific question; expand literature reading, simulation, baselines, or reporting only when the evidence calls for it. Preserve effective work. Do not convert a research process into a mandatory checklist or silently broaden the task.

## Load information only when needed

- **Research standards:** `project_records/long_term_requirements.md` contains the full approved goals, reasoning approach and quality criteria. Consult its relevant sections when doing substantive research design, evaluation or resolving a scientific-standard question; do not reload the entire file for every prompt. Its detailed process is guidance for judgment, not an obligatory sequence of deliverables.
- **Stage-specific guidance:** When the task enters a particular research stage or needs that stage's criteria, read only its matching file under `project_records/research_workflow/`: `01_literature.md`, `02_problem.md`, `03_method.md`, `04_experiments.md` or `05_writing.md`. Do not pre-read the other stages.
- **Project handoff:** Read `project_records/HANDOFF.md` when taking over in a new conversation and the previous research state needs restoring, or when the user explicitly requests it. Do not reread it during an ongoing conversation with sufficient context. Historical entries, including old "next steps", are not current instructions; prefer the latest applicable evidence and the current user request.
- **Project evidence:** Open particular `research/`, `experiments/`, `literature/`, decisions, or Git-history records only to answer the task at hand. Do not sweep the whole repository by default.
- **Tools:** Read `project_records/tool_rules/codex_cli.md` only when invoking or delegating to Codex CLI; read `project_records/tool_rules/git_sync.md` when preparing repository modifications or a completed-checkpoint sync. No tool-use rule should trigger unrelated research.

## Project safeguards and communication

- Preserve `00_基线快照_2026-10-14/` byte-for-byte; preserve scientific meaning, original evidence, equations, assumptions, references and raw data unless the user explicitly directs a justified change.
- Record scientific conclusions in the relevant research or experiment document, key decisions briefly in `project_records/decision_log.md`, and Codex/ARIS prompts and final outputs in `project_records/agent_runs/` when used. Git history records changes; do not create redundant history files.
- Commit and push completed project changes under the standing repository authorization, following the on-demand Git sync rule. Do not force-push.
- Report the concrete question, actual work, results and limits clearly and concisely. On Windows use `python -X utf8`. Do not translate existing material unless requested.
