# AGENTS

This guide exists for agents working in the complex periodic variable repository.

## Context

- This repository is an amalgam of scripts and pipelines used to study a specific class of variable star, complex periodic variables.

## Coding standards

- Write Python that passes PEP 8; keep lines readable and use descriptive names.
- Use Google-style docstrings for every public module, class, function, and method to describe arguments, returns, raises, and any non-obvious behavior.
- Use f-strings for logging messages.
- Favor small, testable functions and keep high-level routines thin when possible.

## Deployment

- During testing and development, this code is deployed on local laptops and servers.  IMPORTANT: during development, use the `cpv` conda environment for all Python development, tests, and automation runs.  For instance, use /Users/luke/local/miniconda3/envs/cpv/bin/python and miniconda3/envs/cpv/bin/pytest directly.  If it doesn't exist, make a suitable uv replacement, also called "cpv".

## Workflow expectations

- Keep `TODO.txt` as the authoritative roadmap: update it whenever new tasks are identified or completed, and mirror finished work in `DONE.txt` to maintain a running history.
- Do not rewrite history (no rebases on published work) unless explicitly requested.

## Testing and quality gates

- Use `pytest` for automated testing; add or update tests whenever behavior changes.
- Before committing, run the relevant subset of `pytest` locally inside the `cpv` environment and fix regressions.

## Additional guidance

- Prefer readable logging and error handling over silent failures; this service has multiple daily stages where diagnostics matter.
- Document assumptions in code comments only when the implementation is non-obvious or relies on domain details (avoid noise).
- Coordinate deployment steps through GitHub pushes; no extra tooling is currently defined.
