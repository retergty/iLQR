---
name: git-update-summary
description: Generates Chinese update summaries by comparing a user-provided git commit id with the current HEAD commit. Use when the user asks to compare a commit with HEAD, summarize changes between commits, generate release notes, list new features, bug fixes, or code line additions/deletions.
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
   - Use the diff content as the source of truth.
   - Use commit messages only as supporting context.
   - For large diffs, inspect the highest-impact files first from `--stat`, then read focused diffs by path.
7. Classify changes into:
   - New features or capabilities.
   - Bug fixes or corrected behavior.
   - Refactors, tests, documentation, or build/tooling changes when relevant.
8. Calculate total code line changes from `--numstat`:
   - Sum numeric added and deleted columns.
   - Treat binary files marked with `-` as binary changes and report them separately instead of counting them as lines.

## Output Format

Default to Chinese unless the user asks otherwise. Use this structure:

```text
从 <commit> 到 HEAD 的更新信息：

提交范围：
- <short-hash> <commit subject>
- ...

新增功能/特性：
- <基于代码 diff 总结的功能点>

修复的 BUG：
- <基于代码 diff 和提交信息总结的修复点>
- 如果没有明确修复，写“未发现明确的 BUG 修复项。”

其他变化：
- <重构、测试、文档、构建配置等>

代码行数变化：
- 新增：<added> 行
- 删除：<deleted> 行
- 净变化：<added - deleted> 行
- 二进制文件：<如有，列出数量或文件名>
```

## Summary Guidelines

- Prefer concise bullet points, each describing user-visible behavior or important internal capability.
- Do not invent features, fixes, tests, or performance improvements that are not visible in the diff or commit messages.
- Mention affected modules when useful, such as `OptimalControl/`, `Approximation/`, `AugmentedLagrangian/`, `TestTools/`, or `Documents/`.
- Distinguish a true bug fix from a refactor or cleanup. Use “修复” only when the diff corrects wrong behavior, build failure, invalid math, unsafe state, or documentation error.
- If the provided commit is not an ancestor of `HEAD`, still compare trees with `git diff <commit>..HEAD`, but state that the revisions are not in a direct ancestor relationship and that the commit list uses `<commit>..HEAD` reachability semantics.
- Do not run `git commit`, change files, or modify git history while generating the update summary.
