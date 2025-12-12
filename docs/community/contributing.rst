Contributing
==================

We welcome contributions from the community to help improve GeodePy! Whether you're fixing bugs, 
adding new features, or enhancing documentation, your input is valuable. Open source projects 
live and die based on the support they recieve.

This document outlines some of the guidlines and advice for contributing to GeodePy.

Code of Conduct
----------------

By participating in this project, you agree to abide by the 
`Python Software Foundation Code of Conduct <https://policies.python.org/python.org/code-of-conduct/>`_. 
Please read it to understand the expectations for behavior when contributing to this project.

Coding Style Guide
------------------

GeodePy uses `Black <https://github.com/psf/black>`_ to keep coding style consistent while still being accessable to all. 
Black uses `PEP 8 <https://peps.python.org/pep-0008/>`_ coding style, an industry standard for python code. Before 
any commits to GeodePy ensure Black has been used.

.. _code:

Code Contributions
------------------

When contributing code please follow these steps:

1. Fork the repository on `GitHub <https://github.com/GeoscienceAustralia/GeodePy>`_.
2. Run tests on current code to ensure it works on your system (See :ref:`Testing <testing>`)
3. Create tests that demonstrate your bug or feature.
4. Make changes, ensuring coding sytle guide is abided by.
5. Run all tests again including one added and ensure all tests pass.
6. Send a Github Pull Request to the repository's **master** branch

Our project maintainers have the last word on if contributions are suitable or not. If your contribution is rejected dont despair!
Following the guidlines above will give you the best chance of getting accpeted.

Documentation Contributions
---------------------------

Documentation imporvements are always welcome! We understand that good documentation is important for all users of a package.
The documentation files can be found in the docs/ folder. They are written in `reStructedText <http://docutils.sourceforge.net/rst.html>`_, 
and use `Sphinx <http://sphinx-doc.org/index.html>`_ to generate the documentation.

When contributing documentation please follow the style of current documentation, having a semi-formal yet friendly approach. 
Ensure any code in documentation is well commeneted to ensure parameters are well understood.

Bug Reports
------------

We welcome all bug reports! Before you raise one though please check the `GitHub issues <https://github.com/GeoscienceAustralia/GeodePy/issues>`_, 
both open and closed, to confirm the bug hastn been reported before. If you do submit a bug report ensure that the bug is clearly described, 
giving the situation that caused the bug and some repeatable code for testing.

Feature Requests
----------------

If you believe a feature is missing, feel free to raise a feature request. Keep in mind that being an open source project requested features may 
or may not be implemented. If there is a feature you really need consider creating it yourself and :ref:`submitting the code <code>`.