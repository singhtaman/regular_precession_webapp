# CLAUDE.md — Physics App Modernization & Vercel Migration

## 1. Mission

Your primary goal is to identify the existing physics application in this repository that is currently deployed on Render, understand how it works, and rebuild or refactor it into a polished, maintainable, production-ready application deployed on Vercel.

The new application must provide a substantially better user interface, better code organization, and an improved user experience while preserving the original application's physics, mathematical models, numerical methods, and scientific correctness.

A reference Vercel application development guide is provided in `VERCEL_APP_GUIDE.md`. Read that guide before making architectural or UI decisions and use it as a reference for building the new application.

**The reference guide defines how to build a high-quality Vercel application. The existing physics application defines what the application must calculate and how its physics must behave.**

Do not confuse these responsibilities.

## 2. Mandatory Initial Investigation

Before modifying any code, inspect the repository thoroughly.

### Identify the existing application

- Inspect the repository structure, README files, configuration files, dependency manifests, and existing documentation.
- Identify the application that is currently deployed on Render.
- Look for Render configuration files, deployment references, service definitions, environment variables, and existing application URLs.
- Determine the framework, programming language, frontend, backend, APIs, scientific libraries, and deployment architecture.
- Identify the application entry points, primary pages, user workflows, computational modules, and data flow.
- Determine whether the application is a monolith, a frontend/backend combination, or a collection of services.
- Identify the existing tests, scientific validation scripts, and example inputs and outputs.

Do not assume the Render application is located at the repository root. Search subdirectories and related packages when necessary.

If multiple applications or deployment configurations exist, identify the correct physics application before proceeding.

### Understand the physics implementation

Locate and document:

- Governing equations and mathematical models.
- Physical parameters, constants, units, and conventions.
- Numerical integration methods, solvers, optimization procedures, and algorithms.
- Initial conditions, boundary conditions, and assumptions.
- Coordinate systems, transformations, and sign conventions.
- Input validation, parameter constraints, and error handling.
- Scientific calculations, intermediate results, and outputs.
- Plotting, data visualization, and result interpretation.
- Any external scientific packages or compiled computational dependencies.

Trace the complete path from user input to numerical calculation and final output.

Create an internal understanding of which code implements the physics and which code implements presentation, routing, infrastructure, or deployment.

Do not begin a large rewrite until you understand the existing application's behavior.

## 3. Reference Vercel App Development Guide

Read `VERCEL_APP_GUIDE.md` in full before implementing the new application.

Use it as the primary reference for:

- Application architecture and project organization.
- Frontend component design and reusable UI patterns.
- Styling, typography, spacing, color systems, and responsive layouts.
- Navigation, forms, controls, and user interactions.
- Charts, scientific visualizations, and results presentation.
- Loading states, error handling, and accessibility.
- API design and communication between frontend and backend.
- Code quality, testing, maintainability, and deployment conventions.
- Vercel configuration and production-readiness practices.

Inspect the reference repository or additional reference files when they are available and authorized for access.

Do not blindly copy the reference application's domain-specific functionality. Reuse its relevant engineering and design patterns while adapting them to the physics application's requirements.

If the reference guide conflicts with the existing physics implementation, preserve scientific correctness and adapt the application architecture accordingly.

If the guide is unavailable, explicitly identify that limitation and proceed cautiously using the existing repository's evidence. Do not pretend to have inspected files you cannot access.

## 4. Non-Negotiable Requirement: Preserve the Physics

**Improve the software around the physics without unintentionally changing the physics itself.**

The existing application is the scientific source of truth.

Unless explicitly authorized by the user, do not change:

- Governing equations or physical models.
- Physical constants, parameter definitions, or units.
- Numerical algorithms, solver settings, tolerances, or convergence criteria.
- Initial conditions, boundary conditions, or assumptions.
- Coordinate transformations or sign conventions.
- Sampling resolutions or discretization choices.
- Formula implementations, computational ordering, or precision in ways that change results.
- The interpretation or meaning of existing outputs.

You may refactor code for clarity, modularity, testability, and maintainability, but any change that could affect numerical results must be investigated and validated.

Prefer extracting existing scientific calculations into well-defined modules or services rather than reimplementing them from scratch.

Do not replace a scientific algorithm with an approximate alternative simply because it is easier to implement in the new frontend or deployment environment.

Do not hardcode expected results, fabricate data, or substitute mock calculations for the real physics.

If the existing implementation contains apparent bugs, questionable equations, inconsistent units, or numerical instability, document the issue and explain its impact. Do not silently alter the scientific model. Seek clarification before making a scientifically meaningful change.

## 5. Establish a Baseline Before Refactoring

Before changing the application, establish a baseline for its current behavior.

Where feasible:

1. Identify representative valid, boundary, and invalid inputs.
2. Run the existing application or computational modules.
3. Record numerical outputs, units, array dimensions, and relevant intermediate values.
4. Save reproducible input/output fixtures or reference results.
5. Identify existing visualizations and user workflows.
6. Record any known limitations or existing failures.

Create regression tests around the scientific calculations whenever practical.

For floating-point calculations, use appropriate absolute and relative tolerances based on the existing numerical methods and scientific requirements. Do not require bitwise equality where it is not scientifically or numerically appropriate.

Compare the new application against the baseline using identical inputs and equivalent execution settings.

A successful UI redesign must not be treated as proof that the physics is correct.

## 6. Redesign the User Interface

Build a significantly more polished and intuitive application, guided by `VERCEL_APP_GUIDE.md`.

The redesign should prioritize:

- A clear information hierarchy.
- A modern, consistent visual design.
- Intuitive navigation and user workflows.
- Well-designed forms and parameter controls.
- Clear labeling of physical quantities and units.
- Useful validation messages and actionable errors.
- Readable equations, numerical results, and scientific notation.
- High-quality, interactive scientific plots where appropriate.
- Responsive layouts for desktop, tablet, and mobile.
- Accessible controls, sufficient contrast, and keyboard navigation.
- Clear loading, processing, empty, and failure states.
- Appropriate precision and significant figures in displayed results.

Organize the interface around the actual scientific workflow. Users should be able to understand what parameters to provide, run a calculation, interpret the results, and inspect relevant plots without unnecessary friction.

Preserve important existing functionality even if its presentation changes.

Do not remove advanced controls, outputs, or scientific features simply to make the interface look simpler.

Avoid superficial redesigns that only change colors and spacing. Improve information architecture, component reuse, interaction design, and the overall usability of the application.

## 7. Improve Code Quality

Use the existing codebase as the starting point and make deliberate, evidence-based architectural improvements.

- Separate scientific computation from presentation and deployment concerns.
- Keep numerical calculations deterministic wherever possible.
- Create reusable UI components and clearly defined module boundaries.
- Use descriptive names, explicit interfaces, and appropriate type checking.
- Avoid duplicated formulas and duplicated business logic.
- Centralize validation and error handling where appropriate.
- Keep configuration and secrets out of source code.
- Remove genuinely obsolete code only after confirming it is unused.
- Preserve important dependencies and computational capabilities.
- Add tests for critical calculations, APIs, and user workflows.
- Document complex scientific algorithms and non-obvious assumptions.
- Avoid unnecessary frameworks, dependencies, and architectural complexity.

Do not rewrite the entire codebase simply to use a preferred framework or language. First determine what can be reused, what should be refactored, and what genuinely needs to be replaced.

## 8. Migrate the Application to Vercel

The target deployment platform is Vercel.

First determine whether the current application can run directly on Vercel or requires a separate computational service.

Inspect the existing Render deployment architecture and identify all functionality that depends on its runtime, filesystem, background processes, persistent storage, or server capabilities.

Choose an architecture that preserves the full functionality of the application.

Depending on the existing implementation, this may involve:

- Deploying a compatible frontend on Vercel.
- Migrating compatible backend endpoints to Vercel Functions.
- Adapting routing and API communication.
- Moving persistent data to an appropriate external database or storage service when required.
- Retaining or replacing computational services that cannot run within Vercel's runtime constraints.
- Moving environment-specific settings into Vercel environment variables.
- Configuring build commands, output directories, runtime settings, and deployment configuration.
- Configuring domains, CORS, and authentication where applicable.

Do not assume that all Render services can be migrated directly to Vercel.

Pay particular attention to scientific workloads that may require substantial CPU time, memory, native libraries, long-running computations, background workers, or persistent processes. Check the current platform limits and compatibility requirements before choosing an implementation.

If the physics calculations cannot run reliably within Vercel's execution limits, preserve the computational engine in an appropriate separate service and deploy the user-facing application to Vercel.

Do not sacrifice scientific correctness or application reliability just to force every component onto Vercel.

Never commit API keys, credentials, or secrets. Do not invent production environment variables or claim that a deployment is configured when it has not been verified.

## 9. Implementation Workflow

Follow this sequence.

### Phase 1 — Discover

- Inspect the repository and locate the Render application.
- Read `VERCEL_APP_GUIDE.md`.
- Map the existing architecture, physics engine, UI, dependencies, and deployment setup.
- Identify constraints and risks.
- Produce a concise implementation plan.

### Phase 2 — Baseline

- Run the existing application and tests where possible.
- Capture representative scientific outputs.
- Identify the essential features and user workflows that must be preserved.
- Establish regression tests for important calculations.

### Phase 3 — Architecture

- Decide what to reuse, refactor, replace, or separate.
- Define the target Vercel deployment architecture.
- Confirm that scientific dependencies and computational workloads are compatible.
- Implement the smallest reliable architectural changes needed.

### Phase 4 — UI and Application Development

- Implement the redesigned interface.
- Refactor components and improve code organization.
- Integrate the actual scientific engine.
- Preserve existing capabilities and input/output semantics.
- Add validation, loading states, and error handling.

### Phase 5 — Verification

- Run available automated tests.
- Run physics regression tests against the baseline.
- Verify representative inputs, edge cases, units, and numerical outputs.
- Verify all major user workflows and visualizations.
- Run the production build.
- Resolve build failures, broken imports, runtime errors, and configuration problems.
- Verify API and backend integration in the target deployment environment where possible.

### Phase 6 — Deployment

- Configure Vercel for the selected architecture.
- Verify environment variables and runtime compatibility.
- Deploy using the available authorized workflow.
- Test the deployed application, including its calculations and critical user flows.
- Document any remaining limitations or steps requiring user action.

Do not claim that a build, test, or deployment succeeded unless you actually ran it and verified the result.

## 10. Working Rules

- Inspect before editing.
- Plan before large refactors.
- Preserve existing scientific behavior by default.
- Use the reference guide for software architecture and UI quality.
- Prefer incremental, reviewable changes.
- Reuse existing dependencies and code when appropriate.
- Avoid unrelated modifications.
- Do not delete important files without understanding their purpose.
- Do not fabricate credentials, test results, scientific data, or deployment status.
- If a necessary decision is ambiguous and materially affects scientific correctness or the deployment architecture, explain the trade-off and ask before proceeding.
- If a minor implementation detail is ambiguous, make a reasonable decision based on the repository and document it.

## 11. Definition of Done

The project is complete when:

- The original Render application has been identified and its behavior understood.
- The reference Vercel guide has been applied where relevant.
- The new UI is materially better and supports the important existing workflows.
- The codebase is more maintainable and appropriately modular.
- The original physics and numerical behavior have been preserved and validated against a baseline.
- Important existing functionality has not been silently removed.
- Tests and production builds pass, or remaining failures are explicitly documented.
- The Vercel deployment architecture is compatible with the application's actual computational requirements.
- The application is deployed and verified on Vercel, or any deployment steps that require additional access or user action are clearly identified.

## Final Principle

**Build a better application, not a different scientific model.**

Use the reference repository to learn how to build and deploy the software. Use the existing Render application to understand and preserve the science. Deliver a polished, reliable Vercel application that combines both.