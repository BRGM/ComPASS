
Coding/Versioning conventions
=============================

Do not hesitate to start discussions on `GitLab <https://gitlab.inria.fr/compass/v4/ComPASS>`_ using issues or merge request. This should be the preferred way.

Versioning conventions
----------------------

Use WIP (Work In Progress) in merge request if you want to prevent an effective merge.

Commits should be carefully formatted (from `Chris Beams <https://chris.beams.io/posts/git-commit>`_):


* Separate subject from body with a blank line
* Limit the subject line to 50 characters (be clear and concise: straight to the point)
* Capitalize the subject line
* Do not end the subject line with a period
* Use the imperative mood in the subject line
* Wrap the body at 72 characters
* Use the body to explain what and why vs. how

The first word of the commit subject may give hints about the nature of the commit (\ *e.g.*\ ):


* WIP: Work In Progress
* BugFix: Commit corresponding to a bug fix. If any try to reference the issue corresponding to the bug.

You can add references in issue or merge request descriptions or in comments. This will update the issue with info about anything related :


* to reference an issue: #123
* to reference a merge request: !123
* to reference a code snippet: $123

Note that issues can be automatically closed through commits using a sepcial syntax (cf. `gitlab documentaiton <https://docs.gitlab.com/ee/user/project/issues/automatic_issue_closing.html>`_\ ).

Branch names should begin by the name of their owner. Don't be shy on branch naming as branch are to be deleted when merge request to the default branch is closed.

Inside the code
---------------

Flags which are pervasive throughout the code and can be located using a grep utility :


*
  `FIXME`: flags a known (potential) bug : None of these should remain in release versions!

*
  `CHECKME`: flags an uncertain development (send a WIP merge request and start a discussion)

*
  `WARNING`: flags a tricky point at code level (ideally those warnings should also be reflected in code documentation)

*
  `OPTIMIZE`: flags a possible enhancement

All flags can be combine with question mark whenever one is not sure about the status of his assertion.

Code specific conventions
-------------------------

When possible, use the prefix ``tmp_`` for the files or directories you don't want to track. They won't be reported by ``git status`` since ``.gitignore`` contains the filter ``tmp_*``.

Fortran
^^^^^^^


*
  prefer lower case fortran

*
  use `fprettify <https://github.com/pseewald/fprettify>`_ for code formating with 3 spaces indentation
