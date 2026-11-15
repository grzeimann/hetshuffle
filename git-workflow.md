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

## Git configuration

Configure your author information once on each computer. Add the following to
your global `~/.gitconfig`, or run the equivalent `git config --global`
commands:

```ini
[user]
    name = Greg Zeimann
    email = gregz@astro.as.utexas.edu

[core]
    editor = nano
```

Check the effective configuration with:

```bash
git config --global --list
```

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

Use one of these commands for a new local checkout:

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

### Commit and publish the change

Run the project's relevant tests or checks before committing. Then review the
staged diff and create a focused commit:

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
