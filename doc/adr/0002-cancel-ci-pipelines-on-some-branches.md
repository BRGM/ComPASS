# 2. Cancel CI pipelines on some branches

Date: 2022-01-11

## Status

Accepted

## Context

We want to be able to skip CI pipelines on some branches to spare limited gitlab runner resources.

Skipped pipeline must not appear as successful, otherwise we could accidentaly merge some modifications
that have not been tested. The remaining options are that the pipeline is cancelled or fails.

Pipelines can be cancelled manually (potentially tedious) or using the gitlab API
(e.g. https://gitlab.com/gitlab-org/gitlab/-/issues/292816). To do so, either the project must be public so that we can
use the CI_JOB_TOKEN variable, either we must expose a TOKEN as a project variable. Yet, this token would potentially
grant any people triggering pipeline rights on the project that they are not supposed to have, with an increased
risk of accidental data loss. Setting the project public would be more satisfying here.

It's quite easy to make pipelines fail depending on the result of a test on the branch name.
The major drawbacks are:
- in the gitlab pipeline landing page we cannot distingush pipelines that fail intentionnaly from other pipelines,
- depending on its own gitlab settings the developer pushing the pipeline will receive a lot of "pipeline failure" notifications.

## Decision

We make the pipeline fail in the first .gitlab-ci.yml job if the branch name starts with "NOCI" (upper case).

## Consequences

Depending on its own gitlab settings the developer pushing the pipeline will receive a lot of "pipeline failure" notifications.
