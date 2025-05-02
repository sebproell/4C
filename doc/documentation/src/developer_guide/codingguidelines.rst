.. _coding-guidelines:

Coding Guidelines
==================

General Guidelines regarding Coding in C++
--------------------------------------------

**Preamble:** |FOURC| has a long history. Large parts were written pre-C++11 and are not modern C++.
In addition, code review was not widely practiced back in these days.
Thus, it is recommended to critically examine the code and design while working on the code base and refer to
general guidelines on C++ coding.
We like to follow the `C++ Core Guidelines <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`_.

The current C++ standard in |FOURC| is C++20 and use of C++20 features is encouraged.

|FOURC|-specific coding guidelines
------------------------------------

In addition, we have a few guidelines that are specific to |FOURC| and are often a product of its long history.
Some of these guidelines are partly covered by the C++ Core Guidelines but are repeated here to emphasize their
importance for the future development of |FOURC|. In these cases, we want to make clear that we are aware that the
code base is not perfect, that we agree with the guidelines, and that we believe the guidelines should be followed
in the future.

**New code should always follow the following guidelines. Touching old code is always a great opportunity to improve it with respect to these guidelines.**

.. _avoid-define-flags:

Avoid define flags
~~~~~~~~~~~~~~~~~~

Do not introduce new define flags. They complicate testing and lead to untested code.

- If you have additional debug code, put it behind ``FOUR_C_ENABLE_ASSERTIONS``.
- If you want to selectively print information during execution of an algorithm, make this a configurable parameter of
  the class or function.

**Exception:** A few flags are defined at configuration time to tell which features are enabled or disabled. These
flags are all managed via CMake.

Avoid header-in-header inclusion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Avoid and actively resolve header-in-header inclusion to speed up compilation time. Use forward declarations instead.

**Attention:** Do not forward declare classes from external libraries such as Trilinos all over the code base. If you
want to avoid including an (expensive) header file via a forward declaration, put that forward declaration in a header
file and include this header. This way, the forward declaration is only in one place.

Use of smart pointers
~~~~~~~~~~~~~~~~~~~~~

If necessary, use smart pointers for their memory management capabilities, e.g., ``std::shared_ptr``, ``std::unique_ptr``
(you may find more information on passing smart pointers
`here <https://www.modernescpp.com/index.php/c-core-guidelines-passing-smart-pointer/>`_

By default, pass parameters by const reference and only expose smart pointers if memory management is necessary.

See also `<https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rr-owner>`_.

Use of ``Teuchos::ParameterList``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prefer parameter container classes over ``Teuchos::ParameterList`` for passing parameters to functions. Often,
a simple struct with all necessary data as public, non-const members is sufficient.

Reason: ``Teuchos::ParameterList`` is a flexible container for input(!) parameters, but it is also very slow and heavy.
In addition, the string-based lookup is not compile-time checked.

Const-correctness
~~~~~~~~~~~~~~~~~

Make code const-correct. Const-correctness is mainly relevant for function parameters and member functions.

**Note:** The member fields of a struct or class should often not be ``const``, as this would prevent the struct or
class from being moved or copied. Instead, use ``const`` member functions to access the fields in a ``const`` context.

Enums
~~~~~

- By default, use scoped ``enum class`` instead of unscoped ``enum``.
- Use enums instead of strings for parameters that take one from a fixed set of values.
- You do not have to write conversion functions between enum constants and strings.
  Use ``EnumTools::enum_name`` to get the name of an enum constant and ``EnumTools::enum_cast`` to convert a string to
  an enum constant. ``EnumTools::operator<<`` and ``EnumTools::operator>>`` allow for easy interaction
  with std streams.

See also `Enum.3 <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#enum3-prefer-class-enums-over-plain-enums>`_.


|FOURC|-specific Naming Conventions
-----------------------------------

- **Namespaces** use CamelCase: ``LinAlg::``
- **class**, **struct** and **enum** types use CamelCase: ``ClassName``, ``StructName``, ``EnumName``
- **Functions** use snake_case: ``function_with_descriptive_name()``
- **Variables** and **enum values** use

    - snake_case, i.e. all lower-case letters separated by underscores, e.g. ``variable_with_descriptive_name``
    - camelCase starting with a lower-case letter, e.g. ``variableWithDescriptiveName``

- **Private class members** end with an underscore: ``variable_``
- **Define flags** are all CAPS: ``FOUR_C_ENABLE_ASSERTIONS`` (please :ref:`avoid define flags<avoid-define-flags>`)

The full specification of the naming convention is listed in the
`.clang-tidy configuration file <https://github.com/4C-multiphysics/4C/blob/main/.clang-tidy>`_.

Variable names must not be just a single letter, because they are impossible to find in a global search operation.
(Exception: loop indices such as i, j, but remember that even loop indices could/should have descriptive names.)


Directory structure
--------------------
The |FOURC| code comes with documentation, example input files and
support scripts. The important subdirectories are the following:

.. list-table::
   :header-rows: 1

   * - Directory
     - Description
   * - ``apps``
     - Contains executables built on top of the |FOURC| library code.
   * - ``cmake``
     - Contains the CMake configuration files to build |FOURC|.
   * - ``dependencies``
     - Contains installation scripts for the external dependencies of |FOURC|.
   * - ``doc``
     - Contains the sources of the documentation you are currently reading.
   * - ``docker``
     - Contains the setup files for docker images for running |FOURC| inside of it.
   * - ``presets``
     - Contains preset files for CMake.
   * - ``src``
     - Contains the bulk of the |FOURC| code organized into several modules.
   * - ``tests``
     - Contains end-to-end tests that run |FOURC| as a whole with optional pre- and post-processing steps. The tests are not directly related to one specific feature (these reside next to the sources).
   * - ``tests/input_files``
     - Contains various valid and running input files. These input files are also used for automated testing.
   * - ``unittests``
     - [Legacy] Contains the source files for stand-alone tests aka unittests. These tests will move closer to the source code into the respective modules.
   * - ``utilities``
     - Contains useful scripts to develop and test |FOURC|.

The top-level directory should only contain files that are commonly expected in this location,
such as the README, LICENSE, and configuration files for tools like clang-format and clang-tidy.

Details of the ``src`` directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``src`` directory contains the modules that make up the |FOURC| library. Each module is
split into a ``src`` and ``tests`` directory. The ``src`` directory contains the production code
of the module, while the ``tests`` directory contains the module-related tests.

Types of source files used within |FOURC|
"""""""""""""""""""""""""""""""""""""""""

The code base consists almost exclusively of C++ source files. To better indicate the intended usage
of a file, both to the build system and developers, the following file extensions are used:

.. list-table::
   :header-rows: 1

   * - File Extension
     - Description
   * - ``.cpp``
     - C++ source files. These files are compiled.
   * - ``.hpp``
     - C++ header files. These files are included in other source files.
   * - ``.fwd.hpp``
     - C++ forward declaration header files. These files contain forward declarations of classes and functions that appear over and over again. They are especially useful for external dependencies to isolate the exposed surface area of the dependency. This type of file is usually included in other header files.
   * - ``.templates.hpp``
     - C++ header files which contain additional template definitions. The reason why you might want to separate certain expensive template definitions from the main header file is to speed up compilation times. Sometimes a template definition requires additional expensive includes, which are not necessary for the main header file. In this case, the template definition is best moved to a ``.templates.hpp`` file.
   * - ``.inst.hpp``
     - C++ header files which contain explicit template instantiations. These files are included in the corresponding ``.cpp`` file to instantiate the templates. These files are not strictly necessary as you can also instantiate the templates directly in the ``.cpp`` file. However, by moving the instantiations to a separate file, you can reuse the same instantiations in multiple ``.cpp`` files which contain parts of the implementation of the same classes.
   * - ``.<ext>.in``
     - These files are templates (not in the C++ sense!) requiring additional configuration via CMake. The CMake script will generate the actual file by replacing placeholders in the template file. The generated file will be placed in the build directory with the extension ``.<ext>``.




