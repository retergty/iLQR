---
name: ilqr-code-development
description: Guides code changes in this iLQR project. Use when the user asks to write, modify, refactor, debug, or test code in the iLQR repository.
---

# iLQR Code Development

## Core Rules

When changing code in this repository:

1. Read the relevant existing implementation before editing.
2. Keep changes local to the requested behavior.
3. Preserve the project's static-size design.
4. Do not introduce dynamic memory allocation in core solver paths.
5. Prefer existing module boundaries and naming conventions.
6. Add or update focused tests when behavior changes.

## Project Constraints

This project is an iLQR optimal control implementation derived from OCS2 ideas but rewritten around stricter engineering constraints:

- State dimension, input dimension, horizon length, and constraint dimensions should remain compile-time configuration.
- Prefer fixed-size Eigen types, `std::array`, and preallocated object members.
- Avoid `std::vector`, heap allocation, runtime resizing, `new`, `make_unique`, and `make_shared` in core runtime paths unless the surrounding code already uses them for non-core setup or tests.
- Avoid adding generic framework layers that obscure the iLQR flow.
- Keep constraints modeled through the existing augmented Lagrangian abstractions unless the user explicitly asks for a different design.

## Code Change Workflow

1. Identify the affected module:
   - Solver loop: `iLQR.hpp`
   - Type plumbing: `iLQRDescriptor*.hpp`, `iLQRTypes.hpp`
   - OCP data: `OptimalControl/`
   - LQ approximation: `Approximation/`
   - Dynamics and rollout: `Dynamics/`, `Rollout/`, `Integration/`
   - Constraints and penalties: `Constraint/`, `AugmentedLagrangian/`, `Penalties/`
   - Riccati recursion: `RiccatiEquations/`
   - Tests and reference tools: `Tests/`, `TestTools/`
2. Read nearby tests before changing behavior.
3. Match the local template style and compile-time dimension patterns.
4. Prefer explicit fixed-size aliases from the existing type system.
5. Decide whether the code change affects `Project_Description.md`.
6. After editing, check lints for touched files and run the narrowest relevant tests when practical.

## Testing Guidance

Use focused tests based on the changed area:

- Approximation changes: `ApproximationTest`
- Dynamics or rollout changes: `DynamicsTest`, `LinearSystemRolloutTest`, `testTimeTriggeredRollout`
- Integration changes: tests under `Tests/integration/`
- Riccati changes: `RiccatiTest`
- Penalty or augmented Lagrangian changes: `PenaltyTest`
- Full solver behavior: `iLQREndToEndTest`, `LinearSystemILQRTest`, `OptimalControl/Exp0Test`

Build and test with CMake presets when needed:

```bash
cmake --preset gcc-debug
cmake --build --preset gcc-debug
ctest --preset gcc-debug
```

## Documentation

For project-level context, consult `Project_Description.md`.

For algorithm background, consult:

- `Documents/Iterative_Linear_Quadratic_Regulator.md`
- `Documents/AL_iLQR.md`

When updating documentation, keep the style structured and formal. Do not copy the user's wording directly if they ask for rewritten project documentation.

When code changes alter public usage, module responsibilities, solver workflow, constraints, memory assumptions, build/test commands, or project conventions, update `Project_Description.md` in the same task.

Do not update documentation for purely local bug fixes, formatting-only changes, or tests that do not change behavior or project conventions.
