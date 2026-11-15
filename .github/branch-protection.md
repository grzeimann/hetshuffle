# Recommended `main` branch protection

Branch protection is a GitHub repository setting and cannot be applied by a
commit in this repository. Apply the following rules to `main` in the
repository Settings, using a branch ruleset or the equivalent branch
protection rule:

- Require a pull request before merging.
- Require at least one approving review when practical.
- Dismiss stale approvals when new commits change the reviewed code.
- Require conversations to be resolved before merging.
- Require the `Python 3.12 tests` status check to pass.
- Require the branch to be up to date with `main` before merging.
- Block force pushes and branch deletion.
- Do not require a coverage percentage yet; first establish a stable baseline.

The required status check is the job name in
`.github/workflows/python-tests.yml`. Keep the workflow on `pull_request` and
do not replace it with `pull_request_target`: fork pull requests must execute
without repository secrets.

If repository policy permits trusted maintainers to merge urgent fixes, use an
explicit administrator bypass rather than weakening the normal pull-request
and CI requirements for everyone.
