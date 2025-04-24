.. _doxygen:

Documenting the code with Doxygen
---------------------------------


|FOURC| uses `Doxygen <http://www.doxygen.nl>`__ to generate documentation from annotated source code.
You can find the current documentation `here <https://4c-multiphysics.github.io/4C/doxygen/index.html>`_.

What is Doxygen and why does |FOURC| rely on it?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Doxygen is a documentation generator, i.e. a tool for writing software reference documentation.
The documentation is written within the code, so documentation lives close by to the code that it is documenting.
This makes it relatively easy to keep up-to-date.
Doxygen can cross reference documentation and code, so that the reader of a document can easily refer to the actual code.

In contrast to the high-level documentation you are currently reading, Doxygen generates the more technical
documentation of the code's API and relevant internals.

|FOURC|'s Policy
~~~~~~~~~~~~~~~~

|FOURC| requires Doxygen documentation for all new code.

It is the responsibility of the developer adding the new code to include the relevant documentation.
During code review, reviewers need to check that documentation is present and understandable.

It is also highly recommended and encouraged that any time you come across undocumented entities,
you take the time to document them. When you struggled to understand a piece of code,
you can be sure that the next developer will struggle too, so it can be a great opportunity to improve the code base.
Pull requests that "only" improve documentation of existing code are a welcome contribution!
Our goal is to have a fully documented code base, and the more proactive every developer is about that, the easier it will be to achieve.


How to write useful Doxygen documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

What constitutes useful Doxygen documentation is best explained with an example.

Consider the following poor example of documentation of a constructor:


.. code-block:: cpp

   /**
    * @brief Constructor
    * @param valid_sections list of valid sections
    * @param legacy_section_names legacy sections (not validated)
    * @param comm MPI communicator
    */
   InputFile(std::map<std::string, InputSpec> valid_sections,
             std::vector<std::string> legacy_section_names, MPI_Comm comm);

This version lacks clarity and feels like the author wanted to satisfy a checklist of items to document.
It provides no insight into the purpose or behavior of the function. By contrast, consider the
following improved version:

---

.. code-block:: cpp

   /**
    * Constructs an InputFile that recognizes the specified @p valid_sections. These sections are
    * associated with InputSpec definitions and will be fully validated accordingly.
    *
    * Sections listed in @p legacy_section_names are also accepted for compatibility with older formats,
    * but they are not validated. It is assumed that these sections contain a list of strings.
    *
    * The constructor also supports the following predefined sections:
    * - "INCLUDES": A list of additional files to include.
    * - "TITLE": Free-form content, typically used for metadata or automatically generated information.
    *
    * The MPI communicator @p comm is required, as reading and distributing input is performed
    * collectively across MPI ranks.
    */
   InputFile(std::map<std::string, InputSpec> valid_sections,
             std::vector<std::string> legacy_section_names, MPI_Comm comm);

This version is concise but readable. It conveys the purpose of the function and its parameters
without overwhelming the reader with unnecessary detail. As another even more complete example,
consider the following version:

---

.. code-block:: cpp

   /**
    * @brief Constructs an InputFile for parsing and validating input sections.
    *
    * Initializes the InputFile with a list of known valid sections and optional legacy
    * section names. Valid sections are fully validated using their corresponding InputSpec,
    * while legacy sections are assumed to contain a list of strings and are not validated.
    *
    * @param valid_sections A map of section names to InputSpec objects. These sections are
    * validated based on their specification.
    *
    * @param legacy_section_names A list of accepted section names that bypass validation.
    * Useful for backward compatibility. These sections are treated as simple string lists.
    *
    * @param comm An MPI communicator used to collectively read and distribute input across ranks.
    *
    * In addition to user-defined sections, the following are always accepted:
    * - "INCLUDES": A list of files to include.
    * - "TITLE": Arbitrary content, useful for metadata or script-generated data.
    */
   InputFile(std::map<std::string, InputSpec> valid_sections,
             std::vector<std::string> legacy_section_names, MPI_Comm comm);

This version is very detailed. The interplay of the parameters is described and all parameters are
explained in detail.

The most important takeaways from these examples are

- Write documentation in full English sentences.
- Describe the purpose of the function and its parameters.
- Describe the interplay between the parameters, return value and (if documenting a class member function)
  the internal state of the class.
- It can be useful to include a `@brief` tag to provide a short description.
- Adapt the length of the documentation to the complexity of the function.
  A helper function like ``print(std::ostream&)`` can often be documented in a single sentence,
  while a complex function may require multiple paragraphs and a detailed description of its parameters.

What needs to be documented?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any entity that is declared inside a header (.hpp) file should be documented. An exception are forward
declarations and namespaces declarations, which do not need to be documented.
Note that Doxygen will not generate documentation for entities that are declared in a source file.
However, it can still make sense to document entities declared inside a source file (.cpp) for readers of the source code.


Building the Doxygen Documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Local build using CMake
"""""""""""""""""""""""""""

|FOURC| has defined a custom build target ``doxygen`` (see also our :ref:`list of custom build targets <build4Cwithcustomtargets>`.
In order to create the Doxygen HTML webpage locally, just issue the command

::

    cd <buildDirectory>
    ninja doxygen


This will build the Doxygen documentation in the directory ``<buildDirectory>/doc/doxygen/html/``.
It can be viewed by accessing ``<buildDirectory>/doc/doxygen/html/index.html`` in any browser.


Pipeline build
""""""""""""""""

The Doxygen documentation is also built when you submit a pull request. A merge to the main branch
will trigger the most recent version of the documentation to be built and published
`here <https://4c-multiphysics.github.io/4C/doxygen/index.html>`_


Choosing good names to simplify documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Documenting source code already starts when writing the code (even before writing comments and documentation).
Here are some general remarks and guidelines taken from Robert C. Martin's book
"Clean Code" to facilitate writing easy-to-understand and well-documented code [Martin08]_.

Use Intention-Revealing Names
"""""""""""""""""""""""""""""""

Writing easy-to-understand and well-documented code starts at selecting descriptive and intention-revealing names for software entities.
Ideally, the name of a variable, function, or class should answer all the big questions.
It should tell you why it exists, what it does, and how it is used.

Use Searchable Names
"""""""""""""""""""""""""

Single-letter names and numeric constants have a particular problem in that they are not easy to locate across a body of text.
Global search will deliver more accurate results, if the variable names are unique and searchable.

Use Pronounceable Names
"""""""""""""""""""""""""

Using pronouncable names makes it easier to discuss code with fellow users and developers.
Don't be afraid of using non-abbreviated names!
Typing on a keyboard is quick, and we use a code formatter to take care of formatting and line breaks.

