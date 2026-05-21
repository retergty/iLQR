---
name: git-commit-message
description: Generates Chinese git commit messages for this iLQR repository by reading Project_Description.md and analyzing git diffs. Use when the user asks to write, polish, summarize, or generate a git commit message for current changes.
---

# Git Commit Message

## Workflow

When generating a commit message in this repository:

1. Read `Project_Description.md` first to understand the project architecture, solver workflow, module boundaries, and static-size design constraints.
2. Inspect the current Git state with `git status --short`.
3. Read changes with Git diff commands:
   - Use `git diff --stat` and `git diff` for unstaged changes.
   - Use `git diff --cached --stat` and `git diff --cached` for staged changes.
   - If both staged and unstaged changes exist, distinguish what is staged from what is only in the working tree.
4. Optionally read recent commit style with `git log -5 --oneline` when choosing tone or length.
5. Summarize the intent of the change, not just the changed filenames.

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
3. <修改点三>
4. <修改点四>
```

Use 1 to 4 numbered points. Keep each point concise and behavior-focused.

## Message Guidelines

- Use verbs that match the change: `新增`, `更新`, `修复`, `重构`, `完善`, `移除`, `整理`.
- Mention the most important affected module or concept, such as `约束布局`, `增广拉格朗日`, `Riccati 递推`, `线性二次近似`, `rollout`, or `项目文档`.
- For documentation-only changes, say `新增/完善文档` instead of implying runtime behavior changed.
- For mixed code and documentation changes, make the title cover both only if both are meaningful.
- Do not invent tests, behavior, or modules that are not visible in the diff.
- Do not run `git commit` unless the user explicitly asks to create the commit.
