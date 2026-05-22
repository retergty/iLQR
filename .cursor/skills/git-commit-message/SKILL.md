---
name: git-commit-message
description: Generates Chinese git commit messages for this iLQR repository by reading Project_Description.md, .gitignore, and staged git diffs. Use when the user asks to write, polish, summarize, or generate a git commit message for staged changes.
---

# Git Commit Message

## Workflow

When generating a commit message in this repository:

1. Read `Project_Description.md` first to understand the project architecture, solver workflow, module boundaries, and static-size design constraints.
2. Read `.gitignore` and use it to avoid analyzing ignored build/output folders or other ignored local artifacts.
3. Inspect the current Git state with `git status --short --untracked-files=all`.
4. Generate the commit message from the staged area only, comparing the index against `HEAD`:
   - Use `git diff --cached --stat`, `git diff --cached --name-status`, and `git diff --cached`.
   - If `git diff --cached` is empty, say there are no staged changes to summarize and ask the user to `git add` the desired files or hunks first.
   - For staged new files, rely on `git diff --cached`; they are included once added to the index.
5. Treat unstaged and untracked files as out of scope for the commit message:
   - Use `git diff --stat` only to notice whether there are extra unstaged changes that are not part of the proposed commit.
   - Use `git ls-files --others --exclude-standard` only to notice untracked files that have not been added yet.
   - Mention these excluded changes briefly only when it helps the user avoid accidentally omitting work.
   - Do not include unstaged or untracked content in the commit message unless the user explicitly asks for a working-tree summary.
6. If a file has both staged and unstaged edits, make the message describe only the staged hunks shown by `git diff --cached`.
7. Optionally read recent commit style with `git log -5 --oneline` when choosing tone or length.
8. Summarize the intent of the staged change, not just the changed filenames.

## Project Context To Preserve

This repository implements an iLQR optimal control solver with augmented Lagrangian constraint handling. Commit messages should reflect the affected architectural area when relevant:

- Type and descriptor plumbing: `iLQRDescriptor*.hpp`, `iLQRTypes.hpp`
- Solver loop and settings: `iLQR.hpp`, `DDPSetting.hpp`, `DDPData.hpp`
- Optimal control problem data: `OptimalControl/`
- Linear-quadratic approximation: `Approximation/`
- Dynamics, rollout, and integration: `Dynamics/`, `Rollout/`, `Integration/`
- Constraints, penalties, and augmented Lagrangian logic: `Constraint/`, `Penalties/`, `AugmentedLagrangian/`
- Riccati recursion and search strategy: `RiccatiEquations/`, `SearchStrategy/`
- Tests and validation tools: `Tests/`, `TestTools/`
- Project or algorithm documentation: `Project_Description.md`, `Documents/`

## Output Format

Default to Chinese unless the user asks otherwise.

Provide a ready-to-use commit message in this structure:

```text
<一句话标题>

1. <修改点一>
2. <修改点二>
3. <按当前修改逻辑继续列出必要修改点>
```

Use as many numbered points as the current change logically needs. Keep each point concise, behavior-focused, and avoid splitting one coherent change into multiple artificial bullets.

## Message Guidelines

- Use verbs that match the change: `新增`, `更新`, `修复`, `重构`, `完善`, `移除`, `整理`.
- Mention the most important affected module or concept, such as `约束布局`, `增广拉格朗日`, `Riccati 递推`, `线性二次近似`, `rollout`, or `项目文档`.
- For documentation-only changes, say `新增/完善文档` instead of implying runtime behavior changed.
- For mixed code and documentation changes, make the title cover both only if both are meaningful.
- Do not invent tests, behavior, or modules that are not visible in the diff.
- Do not run `git commit` unless the user explicitly asks to create the commit.
