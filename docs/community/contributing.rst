Contributing
==================

We welcome contributions from the community to help improve GeodePy! Whether you're fixing bugs,
adding new features, or enhancing documentation, your input is valuable. Open source projects
live and die based on the support they receive.

This document outlines some of the guidelines and advice for contributing to GeodePy.

Code of Conduct
----------------

By participating in this project, you agree to abide by the
`Python Software Foundation Code of Conduct <https://policies.python.org/python.org/code-of-conduct/>`_.
Please read it to understand the expectations for behaviour when contributing to this project.

Coding Style Guide
------------------

GeodePy uses `Ruff <https://docs.astral.sh/ruff/>`_ to keep Python formatting and import ordering consistent. Ruff is
configured in ``pyproject.toml`` and is checked in continuous integration.

Before committing Python changes, run:

.. code:: bash

    ruff format .
    ruff check .
    typos .

.. _code:

Code Contributions
------------------

When contributing code please follow these steps:

1. Fork the repository on `GitHub <https://github.com/GeoscienceAustralia/GeodePy>`_.
2. Run tests on current code to ensure it works on your system (See :ref:`Testing <testing>`)
3. Create tests that demonstrate your bug or feature.
4. Make changes, ensuring the coding style guide is followed.
5. Run all tests again including one added and ensure all tests pass.
6. Send a GitHub Pull Request to the repository's **master** branch

Our project maintainers have the last word on if contributions are suitable or not. If your contribution is rejected don't despair!
Following the guidelines above will give you the best chance of getting accepted.

Documentation Contributions
---------------------------

Documentation improvements are always welcome! We understand that good documentation is important for all users of a package.
The documentation files can be found in the docs/ folder. They are written in `reStructuredText <http://docutils.sourceforge.net/rst.html>`_,
and use `Sphinx <http://sphinx-doc.org/index.html>`_ to generate the documentation.

When contributing documentation please follow the style of current documentation, having a semi-formal yet friendly approach.
Ensure any code in documentation is well commented to ensure parameters are well understood.

Bug Reports
------------

We welcome all bug reports! Before you raise one though please check the `GitHub issues <https://github.com/GeoscienceAustralia/GeodePy/issues>`_,
both open and closed, to confirm the bug hasn't been reported before. If you do submit a bug report ensure that the bug is clearly described,
giving the situation that caused the bug and some repeatable code for testing.

Feature Requests
----------------

If you believe a feature is missing, feel free to raise a feature request. Keep in mind that being an open source project requested features may
or may not be implemented. If there is a feature you really need consider creating it yourself and :ref:`submitting the code <code>`.
