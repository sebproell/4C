# BACI-Testing

BACI is tested via GitLab CI/CD pipelines. Tests can be started by schedulers or whenever something changes in the repository. We test different configurations, for a pipeline to pass it has to succeed for all configurations. This document gives a quick overview over the important steps.

## Terminology

- `job`: A series of shell commands that will be executed on a testing machine.
- `pipeline`: A collection of jobs that will be tested. The green check mark next to a commit means, that the pipeline passed. Some pipelines only contain one job.
- `gitlab-runner`: A daemon on a testing machine that collects jobs from GitLab and runs them. Each runner has tags that are used to filter out the jobs sent to that runner. A detailed documentation can be found [here](https://docs.gitlab.com/runner/configuration/advanced-configuration.html).
- `.gitlab-ci.yml`: The configuration file for GitLab testing. In BACI it is located at `tests/testconfig/.gitlab-ci.yml`. A general introduction to the GitLab testing configuration is given in the [GitLab documentation](https://docs.gitlab.com/ee/ci/yaml/).

## How and when are pipelines started?

Every time something changes in the repository (commit, merge, etc.), GitLab checks for a `.gitlab-ci.yml` file in the repository (at the new state). This file will start ONE pipeline and the jobs in that pipeline (if no jobs are to be run, no pipeline is created). In addition to that pipelines can be triggered over the GitLab web interface and by GitLab schedulers.

## Testing in BACI

### LNM and IMCS

Since Baci is used in different environments and configurations, there are some environment variables defined in `.gitlab-ci.yml` that are passed to the gitlab-runner (for example, the variable `CTEST_CONFIGURE_POSTFIX` is used to pass the option `--useDEAL` to the `do-configure` script).

### Daily tests

Each night a full release and a full debug test are started by a GitLab Scheduler on all configurations. If the release test passes, `doxygen` is build and copied to a folder with constant path (if the documentation is automatically loaded to a server it is better if it is always in the same location, the path is given with a variable and can be different for the different configurations).  

### Minimal tests

Every time something changes in the `master` branch, a pipeline is created that performs the code checks and a minimal test. Only if the code checks are successful, BACI is build and the minimal tests are run. The code check can be performed by any gitlab-runner that picks it up, the build and minimal tests are performed on all configurations.

### User triggered tests

Since we allow merges to `master` only for tested commits, the user has to start a pipeline on the commit when submitting a merge request. This can be done under the GitLab web interface: Goto `CI/CD - Pipelines - Run Pipeline`, then select the branch you want to test (the latest commit on that branch will be tested) and push the button `Create pipeline`. This will first perform the code checks and if they pass a full release build and test on all configurations.

### Output

For each job certain testing output is displayed in GitLab (select job under `CI/CD - Pipelines`). For the minimal tests the full output can be viewed online. The full tests produce too much output for the web interface, only the last 200 lines of each failed test are displayed, as well as a summary of the `ctest` call. In each case, the full output of the `ctest` is written to a `.log` file and copied to a local path defined by the `LOG_FILE_PATH` variable.