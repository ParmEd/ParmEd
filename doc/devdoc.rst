Getting Started
---------------

To get started, you should fork the ParmEd Github repository into your own
Github account.  Push changes there, and submit a pull request (PR) with your
changes (no change is too small for a PR!).

This document outlines some of the general guidelines for working on ParmEd to
try and make sure that the code base is kept clean, organized, manageable, and
reliable (i.e., tested).

Test-driven development
-----------------------

Test-driven development is hard, and often irritating, but I've found that it is
*much* better at avoiding bugs down the road than alternatives.  So when you are
writing code for ParmEd, generate a test case *and write at least the beginning
of the test first*.

Then, before you write any code, *run* the test suite to make sure that it
fails! If the test doesn't fail before you've written any code, then your test
is clearly not being run, or not being run correctly (which means there is no
way that your test could prevent a breakage in the future).

Once you have a failing test case, write code until the test case passes. Then
add more test until the test case fails again. Then update the code until the
test passes again. Rinse-and-repeat until the functionality you are adding is
complete and fully-tested.

Coding Style
------------

For the most part, the coding style here follows the style guide `PEP 8
<https://www.python.org/dev/peps/pep-0008/>`_. Of the guidelines listed there,
the following I consider very important (and will block any PRs):

    - Indentation is done with 4 spaces. No source code file should contain
      tabstops, and all indentation must be consistent.
    - Source code should be ASCII encoded
    - Please make sure comments and docstrings are written in English
    - Use only absolute imports or *explicit* relative imports
    - Always use ``print`` as a function (via the ``print_function`` futures
      import) -- see below for why this is a hard requirement.

Things I would like to see, but are not as important:

    - Try to keep lines limited to 80 characters, but do not sacrifice
      readability
    - Class definitions should be done in ``CamelCase`` naming style. Method,
      attribute, and variable names should be done in
      ``underscore_separated_style``. Constants should be named in ``ALL_CAPS``.

Target Python Version
---------------------

ParmEd officially supports Python 2.7 and Python 3.3+ in the same code base
through the use of the ``six`` module (although some earlier versions of Python
*may* work, no effort is made to support Python 2.6 and lower or Python 3.2 and
lower). The ``six`` module is included as a subpackage in the ``parmed.utils``
package, and all components should be imported from there.

In particular, the ``range`` and ``zip`` builtins should be imported from
``parmed.utils.six.moves`` rather than relying on the standard versions. This is
because in Python 3, ``range`` and ``zip`` return efficient iterators, while in
Python 2 they return potentially inefficient lists. The Python 2-equivalent
versions of the ``range`` and ``zip`` iterators are ``xrange`` and
``itertools.izip``, respectively, which are the *actual* functions defined
within ``parmed.utils.six.moves`` for Python 2.

Supporting both Python 2 and Python 3 in the same code base means that you must
never use the ``print`` statement, implicit integer division, or relative
imports (as those will *not* work in Python 3). You are encouraged to put the
following line at the top of each Python module in ParmEd::

    from __future__ import division, absolute_imports, print_function

In addition to supporting CPython 2.7 and 3.3+, ParmEd also strives to support
the ``pypy`` interpreter.  Under the ``pypy`` interpreter, PDB file parsing
is about 3-5x faster *using the same code*. Since ParmEd parses a substantial
amount of metadata from RCSB/wwPDB files, the added efficiency can be useful for
mining the PDB. As such, modules that require CPython (such as ``scipy``,
``scikit-learn``, ``matplotlib``, ... etc.) must be made optional so that the
functionality *not* relying on those modules can be used within ``pypy``.

Tests
-----

Please keep test cases relatively small -- Github likes to keep repositories
under 1 GB in size, so test files should be relatively small. They should also
be kept small enough to run in a short amount of time, but they should also
ensure that the function being tested has as close to full coverage as can be
managed.

The ParmEd tests are run in the ``ParmEd/tests`` directory with the command::

    nosetests -vs .

which will automatically find and run all ``test_*.py`` files. The tests are
implemented using the ``unittest`` module in the Python standard library. Look
at existing tests for examples.

While developing new functionality, I find it helpful to run the tests directly
from the ParmEd directory with the following command::

    nosetests -vs test/test_parmed_structure.py

Obviously in this case, ``test_parmed_structure.py`` is replaced with whichever
test module you are working on. You can select *specific* tests using the ``-m``
flag specifying a regex that matches the test case method.  For example::

    nosetests -vs test/test_parmed_structure.py -m add_atom

will test both the ``test_add_atom`` and ``test_add_atom_to_residue`` methods.
This is an easy way to run tests quickly while working on new methods *without*
having to run ``python setup.py install`` after every change. Note that when you
run tests from the root ParmEd directory, however, the imported ParmEd
repository will not have any Python extensions installed (meaning that the tests
relying on them -- like the test for the Amber optimized reader -- will fail).

ParmEd utilizes the Travis continuous integration server to perform automatic
tests of all pull requests. Tests generally must pass these tests before being
considered for merge into the master branch.

Documentation
-------------

Just as important as the code is *documenting* the code so that users can learn
how to use it. Documentation is tracked as reStructuredText (``.rst``), and
`Sphinx <http://sphinx-doc.org/>`_ is used to render the documentation. Changes
to the documentation can be made part of the pull request adding new
functionality, or it can be its own pull request dedicated entirely to
documentation improvements.

The actual web docs are generated periodically and pushed directly to the
``gh-pages`` branch. This process is not handled by PRs.
