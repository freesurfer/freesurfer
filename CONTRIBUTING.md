# Contributing to FreeSurfer

If you have improvements to the FreeSurfer code base, send us a pull request! If you'd just like to point out a bug, please open up an issue or contact the [help list](https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSupport).

## Developing

Our [developers guide](https://surfer.nmr.mgh.harvard.edu/fswiki/DevelopersGuide) walks through the necessary instructions for correctly configuring and building the FreeSurfer source code. For those new to git, our [GitHub guide](https://surfer.nmr.mgh.harvard.edu/fswiki/GitHub) covers some introductory topics on the suggested workflow.

Commit messages should be descriptive and prefixed by one of the following context tags:

 * `rf:` refactoring
 * `bf:` bug fix
 * `nf:` new feature
 * `doc:` documentation updates
 
So for example, a commit fixing a bug in the surface placement might be titled something along the lines of `bf: fixed vertex sse computation`.

## Submitting changes

For significant changes (anything beyond updating documentation), discuss your proposition via git issues, the help list, or any other method of communication with LCN developers before submitting a pull request. Additionally, make sure you've read our [Code of Conduct](https://github.com/freesurfer/freesurfer/blob/dev/CODE_OF_CONDUCT.md) and run all appropriate unit and regression tests (successfully) as described in the [build guide](https://surfer.nmr.mgh.harvard.edu/fswiki/BuildGuide).

Please document your pull request with details including:

  - An explanation of what you've changed... and *why*
  - References to issues or previous discussions, where appropriate
  - A summary of the tests that you've run on your branch

**Thank you for contributing!**
