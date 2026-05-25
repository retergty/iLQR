---
name: git-update-summary
description: Generates high-level Chinese work reports by comparing a user-provided git commit id with the current HEAD commit. Use when the user asks to compare a commit with HEAD, summarize work content between commits, generate update reports, release notes, feature summaries, bug-fix summaries, or code line additions/deletions.
---

# Git Update Summary

## Workflow

When the user provides a git commit id and asks for update information against the current `HEAD`:

1. Treat the provided commit id as the comparison base and `HEAD` as the target.
2. Verify both revisions before analysis:
   - `git rev-parse --verify <commit>^{commit}`
   - `git rev-parse --verify HEAD^{commit}`
3. Inspect repository state with `git status --short`.
   - The summary must describe committed changes between `<commit>` and `HEAD`.
   - If the working tree has uncommitted changes, mention that they are excluded unless the user explicitly asks to include them.
4. Gather commit messages between the base and `HEAD`:
   - `git log --reverse --format='%h %s' <commit>..HEAD`
   - If this returns no commits but `git diff <commit>..HEAD` is non-empty, explain that the revisions differ by tree content but there are no commits reachable from `HEAD` that are absent from the provided commit.
5. Gather change statistics:
   - `git diff --find-renames --stat <commit>..HEAD`
   - `git diff --find-renames --numstat <commit>..HEAD`
   - `git diff --find-renames --name-status <commit>..HEAD`
6. Read the actual code changes with `git diff --find-renames <commit>..HEAD`.
   - Use the diff content as the source of truth for understanding intent and scope.
   - Use commit messages only as supporting context.
   - For large diffs, inspect the highest-impact files first from `--stat`, then read focused diffs by path.
7. Classify changes into:
   - New features or capabilities.
   - Bug fixes or corrected behavior.
   - Refactors, tests, documentation, or build/tooling changes when relevant.
8. Calculate total code line changes from `--numstat`:
   - Sum numeric added and deleted columns.
   - Treat binary files marked with `-` as binary changes and report them separately instead of counting them as lines.
9. Convert implementation details into report-level language.
   - Explain what work was completed and why it matters.
   - Avoid code snippets, function-level walkthroughs, raw diff details, and overly specific implementation mechanics.
   - Mention modules or subsystems only when they help readers understand the work scope.

## Output Format

Default to Chinese unless the user asks otherwise. Write like a concise work content report, not a code review. Use this structure:

```text
从 <commit> 到 HEAD 的工作内容报告：

概览：
- <根据工作内容多少自动增减要点，概括本区间完成的主要工作方向和整体价值>

新增功能/特性：
- <高层描述新增能力或增强点，不展开具体代码实现>

修复的 BUG：
- <高层描述修复的问题、影响范围或行为改善>
- 如果没有明确修复，写“未发现明确的 BUG 修复项。”

其他变化：
- <重构、测试、文档、构建配置、工程化完善等工作内容>

提交信息参考：
- <short-hash> <commit subject>
- ...

代码行数变化：
- 新增：<added> 行
- 删除：<deleted> 行
- 净变化：<added - deleted> 行
- 二进制文件：<如有，列出数量或文件名>
```

## Summary Guidelines

- Prefer concise bullet points, each describing work outcome, user-visible behavior, module-level capability, or engineering value.
- Let the number of overview bullets follow the amount of work; do not force a fixed count.
- Keep the report at a high level. Do not include code snippets, detailed algorithms, function bodies, line-by-line changes, or patch excerpts.
- Avoid file-by-file narration unless a file represents a major independent deliverable.
- Do not invent features, fixes, tests, or performance improvements that are not visible in the diff or commit messages.
- Mention affected modules when useful, such as `OptimalControl/`, `Approximation/`, `AugmentedLagrangian/`, `Tests/Include/`, or `Documents/`.
- Distinguish a true bug fix from a refactor or cleanup. Use “修复” only when the diff corrects wrong behavior, build failure, invalid math, unsafe state, or documentation error.
- If the provided commit is not an ancestor of `HEAD`, still compare trees with `git diff <commit>..HEAD`, but state that the revisions are not in a direct ancestor relationship and that the commit list uses `<commit>..HEAD` reachability semantics.
- Do not run `git commit`, change files, or modify git history while generating the update summary.
