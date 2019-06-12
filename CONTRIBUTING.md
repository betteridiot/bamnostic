# How to contribute to bamnostic
----
**First and foremost:** by contributing in this project, you agree to abide by
the Contributor Covenant [code of conduct](https://github.com/betteridiot/bamnostic/blob/master/CODE_OF_CONDUCT.md)

## Our Promise:
The maintainers of this project promise to address each issue/question/pull
request in the following manner:
* Prompt acknowledgment of receipt of issue/question/pull request
    * Potentially assigning to a specific maintainer
* If needed, a description of when work on a given issue has started, or an explanation of why the issues/pull request is not being addressed
* Closing the issue upon the maintainer's determination of issue resolution.

## General Questions

* Read the [documentation](https://bamnostic.readthedocs.io/en/latest/?badge=latest) on the project
* Search the GitHub [Issues](https://github.com/betteridiot/bamnostic/issues) for this project to see if this question has already been addressed
* If neither of these avenues answer your question, please feel free to create
a new [issue](https://github.com/betteridiot/bamnostic/issues) that has a 
succinct (but informative) subject line. In the body of the issue, please ask
your question, including any context regarding your specific problem or inquiry.

A *bad* general question would look like:
> "Your module will not work on my computer. How do I fix it?"

A *good* general question would look like:
> "I am running your code on Windows 10 version 1803 in a Anaconda build of Python version 3.6. When I try to \<some implementation\>, I get \<some error\>. Here is an example of how I invoked your code: \<insert code invocation\>. Can you tell what is going on?"

## Bug reporting

* Ensure the bug has not been already reported:
    * Search the GitHub [Issues](https://github.com/betteridiot/bamnostic/issues) for this project
* If the bug has not been previously reported, feel free to open a new issue:
    * Please follow the Bug report template provided on issue creation as closely as possible

## Bugfixes, patches, or documentation corrections
Any additions to the code should follow the [Google Python Style Guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md) or [NumPy Style Guide](https://numpydoc.readthedocs.io/en/latest/) for documentation purposes.

The preferred method of contributing in the form of a pull request is to fork
the latest version of the project (likely the `master` branch), and cloning that
fork to your local machine:

```bash

git clone git@github.com:your-username/bamnostic.git

```

Commit your changes to your fork, and submit a [Pull Request](https://github.com/betteridiot/bamnostic/pulls)

The maintainers will review your PR within a timely manner. Be aware that the 
maintainers may request an improvement or alternative, and also reserve the
right to reject any Pull Request if it does not meet sytle or community
guidelines.

A good PR should contain the following items:
* Follows the [Google Python Style Guide](https://github.com/google/styleguide/blob/gh-pages/pyguide.md) or [NumPy Style Guide](https://numpydoc.readthedocs.io/en/latest/).
* **Tests:** either doctests for unit tests (written using `pytest`)
* Contains a helpful/meaningful **commit message**.
