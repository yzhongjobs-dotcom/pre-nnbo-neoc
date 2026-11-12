# Codex CLI use (read only when invoking Codex)

Claude may delegate concrete implementation, simulation, data analysis or independent checking when useful. Divide work flexibly, but prioritize research quality over model cost or minimizing tokens.

1. Delegate a bounded task: state the exact scientific question, the files/evidence to use, what to implement or verify, and how the result should be checked. Do not pass the entire historical research backlog as an implicit assignment.
2. Choose the available Codex model and reasoning effort according to task difficulty and the need for independent review; do not hard-code a model version for all tasks.
3. Continue an existing Codex session when continuing the same task (`codex exec resume <id>`); use a new session for a genuinely separate task or independent review. Confirm the installed CLI's actual supported commands/options before use.
4. Check Codex's code, mathematical claims and reported results before accepting them. Preserve the raw task prompt and final output in `project_records/agent_runs/`; do not keep process logs by default.
5. Changes to the scientific question, Method architecture or research priorities follow the current user request and shared research principles, not a suggestion made by the delegated agent.
