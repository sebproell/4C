.. _4Ctesting:

Testing
=======

Clean code goes hand in hand with code testing, that is why in |FOURC| we set great value upon the latter.
This page is concerned with everything related to testing in |FOURC|.

Overview on testing mechanisms
------------------------------

|FOURC| uses the ``ctest`` environment for testing (refer to `ctest <https://cmake.org/cmake/help/latest/manual/ctest.1.html>`_)
and relies on a variety of test mechanism to guarantee and maintain its intended behaviour.
Tests fall into these categories:

- Regression tests/end-to-end tests: test a typical simulation run of |FOURC| with an input file
- Unit tests: isolate and exercise specific units of source code independently from other parts

|FOURC| tests can be triggered through various mechanisms:

- On Github, we use Github actions on every pull request and run additional tests every night.
- Locally, one can trigger ctest to run all or some tests.

    - Filter tests by adding ``-R <regex>``. Only tests including ``<regex>`` in their names are executed.
    - Exclude tests by adding ``-E <regex>``. Tests with ``<regex>`` in their names are not executed.
    - Filter tests by adding ``-L <regex>``. Only tests with label ``<regex>`` in their names are executed.
    - Skip the clean up test by adding ``-FC test_cleanup``.

- If you want to run ctest on a single input file inside the ``tests/input_files`` directory,
  run ``ctest -R <input_file>`` where ``<input_file>`` is the input file name.


Guidelines for |FOURC| input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An input file of a CI test (residing in ``<4C_sourcedir>/tests/input_files``) should include:

- A short but informative summary of what is tested, maybe also which results are compared,
  in the header
- No unused sections
- No unused parameters
- No alternative input parameters/values as comment

Executing |FOURC| unit tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Configure and build |FOURC| as described :ref:`here <installation>`.
The |FOURC| unit tests are included in ctest as part of the minimal tests and also in the full test suite:

::

    ctest -L minimal
    ctest -R unittests

The test executables are located in the ``tests/`` directory in the build folder. They support the ``--help`` argument,
which can be used to get a list of available options, e.g. to filter for specific tests.

Many IDEs also come with plugins or support for GoogleTest allowing to run tests directly from the IDE.

How to add unit tests to |FOURC|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|FOURC| uses `GoogleTest <https://github.com/google/googletest>`_ for unit testing. If you are new to this framework,
read the `primer <https://google.github.io/googletest/primer.html>`_.
In general, we recommended to look through existing unit tests directories first to get an idea on how the tests are organized and how GoogleTest can be used.

Unit tests reside close to the module containing the tested functionality, namely in the ``tests`` directory next to a module's ``src`` directory.
It can be a good idea to have one unit test file per source file. Conventionally, this test file is named as the source file with the suffix ``_test``.
Unit tests are picked up automatically by CMake through ``four_c_auto_define_tests()``.