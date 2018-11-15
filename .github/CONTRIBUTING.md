# How to contribute

We welcome contributions from external contributors.
This document describes the contribution process: from proposing a change to
getting it merged into OpenRSP.
The process for contributing is exactly the same for the core development team
and for external contributors:

* Maintainers do not push directly to the repository.
* Maintainers do not review their own patches.
* Approval of one or more maintainers is sufficient for trivial patches.
  Trivial patches are typos and trivial bugfixes
* For any patch that is not trivial, two maintainers need to review and approve the patch.

## Getting Started

* Make sure you have a [GitHub account].
* [Fork] the [openrsp/openrsp] repository on GitHub.
* On your local machine, [clone] your fork of the OpenRSP repository.
* The [Psi4](http://psicode.org/) documentation has more detailed instructions for interacting with your fork which can be found
  [here](http://psicode.org/psi4manual/master/build_obtaining.html#faq-forkpsi4public).
  and [here](http://psicode.org/psi4manual/master/build_obtaining.html#faq-githubworkflow).

## Making Changes

* Add some really awesome code to your local fork. It's usually a [good idea] to
  make changes on a [branch] with the branch name relating to the feature you
  are going to add. A style guide is available in the [STYLE_GUIDE.md file].
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of OpenRSP on GitHub and open a [pull request] (PR)
  __towards the `master` branch__.
  Note that after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.
  Each commit added to the PR will be validated for mergeability, compilation
  and test suite compliance; the results of these tests will be visible on the
  PR page.
* The title of your pull request should be marked with `[WIP]` if it’s a work
  in progress. Feel free to use as many labels as you think necessary!
* Update the [`ChangeLog.rst`] file. We follow the conventions and recommendations detailed at
  [keepachangelog.com]
* When you're ready to be considered for merging, check the "Ready to go" box
  on the PR page to let the OpenRSP team know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration services return passing checkmarks, and maintainers give "Approved" reviews.

## Licensing

* We do not require any formal copyright assignment or contributor license agreement.
* **Any contributions intentionally sent upstream are presumed to be offered under the terms of the OSI-approved [LGPL v2.1 License].**

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
* [Developer Certificate of Origin versus Contributor License Agreements](https://julien.ponge.org/blog/developer-certificate-of-origin-versus-contributor-license-agreements/)


[GitHub account]: https://github.com/signup/free
[Fork]: https://help.github.com/articles/fork-a-repo/
[openrsp/openrsp]: https://github.com/openrsp/openrsp
[clone]: https://help.github.com/articles/cloning-a-repository/
[good idea]: http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/
[branch]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/
[STYLE_GUIDE.md file]: ../STYLE_GUIDE.md
[Developer Certificate of Origin]: https://developercertificate.org/
[pull request]: https://help.github.com/articles/using-pull-requests/
[`ChangeLog.rst`]: ../ChangeLog.rst
[keepachangelog.com]: https://keepachangelog.com/en/1.0.0/
[LGPL v2.1 License]: ../LICENSE
