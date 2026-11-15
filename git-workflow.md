# Git Workflow

This document describes the recommended GitHub workflow for **hetshuffle**.
The goal is to keep `main` stable, make changes easy to review, and make it
straightforward to back out work that is not ready or does not work as
expected.

The proposed testing and CI rollout is documented in
[testing-and-ci-plan.md](testing-and-ci-plan.md).

## GitHub setup

- Sign in to [GitHub](https://github.com/).
- Confirm that you can access the [hetshuffle repository](https://github.com/grzeimann/hetshuffle).
- Configure GitHub authentication if needed. SSH is recommended; see
  [GitHub's SSH key documentation](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

The repository's canonical remote is:

```text
https://github.com/grzeimann/hetshuffle.git
```

For an SSH remote, use:

```text
git@github.com:grzeimann/hetshuffle.git
```

## Contributor access and onboarding

The normal collaboration path is direct work in the canonical repository.
People who are not yet contributors should contact a repository maintainer
before creating a fork or opening a PR. The request should include their
GitHub username and a short description of the work they intend to do.

After confirming the contributor, a maintainer should add them to the
repository with GitHub `Write` access. This allows them to push feature
branches while branch protection continues to prevent direct changes to
`main`. Contributors should not need administrator or maintainer-level access
for ordinary development.

Once the invitation is accepted, the contributor should clone the canonical
repository and use `origin` for the canonical remote. All work still goes on a
feature branch and through a pull request; write access is not permission to
push directly to `main`.

Forks remain a fallback for people who are awaiting access, contributing from
outside the collaboration, or unable to receive repository access. In that
case, use the fork instructions below and open a PR against
`grzeimann/hetshuffle:main`.

## Git configuration

Configure your author information once on each computer. Add the following to
your global `~/.gitconfig`, or run the equivalent `git config --global`
commands:

```ini
[user]
    name = Your Name
    email = you@example.com

[core]
    editor = nano
```

Check the effective configuration with:

```bash
git config --global --list
```

## Create a development environment

hetshuffle supports Python 3.12 or newer. The development extra installs the
runtime dependencies, pytest, coverage, Ruff, and pre-commit:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --editable '.[dev]'
```

Contributors who use Conda can instead create the runtime environment from
`environment.yml`, activate it, and install the editable development extra:

```bash
conda env create --file environment.yml
conda activate hetshuffle
python -m pip install --editable '.[dev]'
```

The standard workflow includes `reproject` and PostgreSQL catalog access.
`pyds9` is optional and is not required for ordinary tests or pull requests.

Install both local hooks after installing the development extra:

```bash
python -m pre_commit install --hook-type pre-commit
python -m pre_commit install --hook-type pre-push
```

The pre-commit hook checks file hygiene and runs Ruff on the test suite. The
pre-push hook runs the complete pytest suite under coverage. CI remains the
authoritative clean-environment check, so do not bypass a failed hook without
understanding and documenting the reason.

Verify the setup before starting work:

```bash
python --version
python -m pre_commit run --all-files
make ci
```

The first pre-commit run may update whitespace or other file-ending details;
review those changes and run it again until it passes.

## Standard feature-branch workflow

1. Start from an up-to-date local `main` branch. Do not develop directly on
   or push directly to `main`.
2. Create one feature branch for each bug fix, feature, or focused change.
3. Commit related changes with clear, imperative commit messages.
4. Push the feature branch and open a pull request (PR) on GitHub.
5. Use the PR discussion for review, testing notes, and follow-up changes.
6. Merge the PR into `main` only after the change is ready and any required
   checks or review are complete.
7. Delete the merged branch locally and on GitHub when it is no longer needed.

### Clone the repository

Approved contributors should use one of these commands for a new local
checkout:

```bash
git clone git@github.com:grzeimann/hetshuffle.git
cd hetshuffle
```

Or, if using HTTPS:

```bash
git clone https://github.com/grzeimann/hetshuffle.git
cd hetshuffle
```

### Start a change

Before creating a branch, make sure any existing local work is committed or
otherwise safely saved:

```bash
git status
git switch main
git pull --ff-only origin main
git switch -c descriptive-branch-name
```

Use short, descriptive branch names such as `fix-config-path` or
`add-target-visualization`. Avoid spaces in branch names.

If you are using the fallback fork workflow, clone your fork and add the
canonical repository as `upstream`:

```bash
git clone git@github.com:YOUR-USER/hetshuffle.git
cd hetshuffle
git remote add upstream https://github.com/grzeimann/hetshuffle.git
git fetch upstream
git switch -c descriptive-branch-name upstream/main
```

For a direct collaborator checkout, use the canonical clone commands above and
start from `origin/main` instead.

### Commit and publish the change

Make the edit, add or update tests when behavior changes, and review the diff.
The pre-commit hook runs automatically when you commit. Run the full local
gate explicitly before publishing:

```bash
make ci
```

Then review the staged diff and create a focused commit:

```bash
git status
git diff
git add path/to/changed-file.py
git diff --cached
git commit -m "Fix configuration path handling"
git push -u origin descriptive-branch-name
```

Push subsequent commits with:

```bash
git push
```

Open a PR from the pushed branch. GitHub Actions will install the package in a
fresh Python 3.12 environment and run the `Python 3.12 tests` check. The PR
template records the exact local commands, scientific impact, and any
data-dependent or manual verification that remains. Keep a new or untested
feature as a draft until those details and the general test results are ready.

For a fork, open the PR from `YOUR-USER:descriptive-branch-name` into
`grzeimann/hetshuffle:main`. Fork pull requests run without repository secrets.

Keep commits understandable and avoid messages such as `did stuff`. If a
change is still in progress, mark the PR as a draft or include `WIP` in its
title. Mark it ready for review once the implementation, tests, and
documentation are in good shape.

### Update a branch with changes from `main`

Before updating a long-lived feature branch, make sure the working tree is
clean. Then fetch the latest remote history and incorporate `main`:

```bash
git status
git fetch origin
git switch descriptive-branch-name
git merge origin/main
```

Resolve any conflicts carefully, run the relevant tests again, and push the
result. If the project adopts a different integration policy later, document
that policy here before using history-rewriting commands such as `git rebase`
or `git push --force`.

## After a PR is merged

Update the local `main` branch and remove the finished feature branch:

```bash
git switch main
git pull --ff-only origin main
git fetch --prune origin
git branch -d descriptive-branch-name
```

The lowercase `-d` deliberately refuses to delete a branch whose changes are
not merged. Use `-D` only after verifying that the branch is no longer needed.

If local work becomes confusing, stop and check `git status` and
`git log --oneline --decorate --graph --all` before taking corrective action.
When appropriate, make a backup branch or fresh clone rather than deleting
unreviewed work.
