# Coding/Versioning conventions

## Versioning conventions

Commits should be carefully formatted :
  * Separate subject from body with a blank line
  * Limit the subject line to 50 characters (be clear and concise: straight to the point)
  * Capitalize the subject line
  * Do not end the subject line with a period
  * Use the imperative mood in the subject line
  * Wrap the body at 72 characters
  * Use the body to explain what and why vs. how
taken from [Chris Beams] (https://chris.beams.io/posts/git-commit/)

The first word of the commit subject may give hints about the nature of the commit:

  * WIP: Work In Progress
  * BugFix

Branch names should begin by the name of their owner.

## Inside the code

Flags which are pervasive throughout the code and can be located using a grep utility :

  * FIXME: flags a known (potential) bug : None of these should remain in release versions!
  
  * CHECKME: 

  * WARNING: flags a tricky point at code level (ideally those warnings should also be reflected in code documentation)
  
  * OPTIMIZE: flags a possible enhancement
  
All flags can be combine with question mark whenever one is not sure about the status of his assertion. **Do not hesitate to start discussions on GitLab using issues or merge request. This should be the preferred way.** Use WIP (Work In Progress) in merge request if you want to prevent an effective merge.





