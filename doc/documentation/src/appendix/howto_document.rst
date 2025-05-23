.. _writingdocumentation:

Writing documentation
=========================

Usage
-----

The documentation is generated by `sphinx <https://www.sphinx-doc.org>`_ from the documentation
files in the ``doc/documentation/src`` directory.

.. figure:: /_assets/documentation_workflow.png
   :width: 500px
   :align: center

   From ascii files to websites

The ascii files are written in restructuredText or markdown.
While restructuredText conversion is the native implementation of sphinx, markdown works well, too.
Both can be converted to an html document and provided to https://www.readthedocs.org.

Here, I'll give some hints, mainly for myself ;-)
You'll find other useful information, e.g., in the `official documentation <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_;
I also like `this doc <https://books.dehlia.in/writing-with-ratatouille/toc/>`_. Otherwise, GIYF.

- You can emphasize text both in markdown and restructuredText by making it **bold** (using double aserisks before and after: ``**bold**``,
  or print it in *italic* (using single asterisk before and after: ``*italic*``).

- **Headers** are sorted automatically. Actually, it is not necessary to use specific header,
  since Sphinx sorts the header by itself. However, it's good practice to use a specific format.
  In general one has to underline the title with the following characters

  - Parts: Both overlining and underlining the title with ``=======``
  - Sections: ``=====``
  - Subsections: ``------``
  - Subsubsections: ``~~~~~~~~``

  .. note::

     The underlines (and the overlines) must be at least as long as the title itself.
     A title (including the underline) has to be surrounded by a blank line.

- **Lists** (as the one addressing the features of reStructuredText here) can be numbered or unsorted.
  Any list must be surrounded by blank lines.
  An unordered list with bullet points is created by dashes like this one::

     Here, we'll start a list

     - Point one
     - Point two
     - Point three

     further text

  A numbered list is created like the unordered ones, simply replace the asterisk with ``# .`` as used in the following.

   #. A list within a list must be indented by three blanks.
   #. Also, before and after the inserted list, blank lines are needed.

- **Links to other headers** (internal links) can easily be included in restructuredText by
  ``:ref:`Descriptive link name <sectionname>``, e.g., :ref:`Workflow <4Cworkflow>`.
  The link target can be a section name, but also an explicit link target.
  The latter is preferred, since there may be duplicate section names in the documentation.
  A link target must be entered by ``.. _targetname:``, i.e., first two dots, then a blank followed by an underline,
  then the name of the link target, and finally a colon.

  From my small experience it seems that the link to a section name does not contain spaces,
  even if the section title does contain them
  (However, I saw spaces in other documentations many times).
  As an example, you'll find a link to the top of this page :ref:`here <writingdocumentation>`.
  Links also work between pages, of course.

- **Links to other websites** are included  by, e.g., ``http://www.readthedocs.org`` (<http://www.readthedocs.org>),
  that is, you may simply type the website and it is recognized as a link.
  If the link should show a different text, use back quotes and an underscore at the end,
  like this: ```ReadTheDocs <http://www.readthedocs.org>`_`` , which then looks like `ReadTheDocs <http://www.readthedocs.org>`_.
  Note that a blank is needed before the opening angle bracket.

- **References** are given in ``references.rst``.
  You'll see how to prepare such a reference as soon as you open the file.
  The corresponding citation is created at any place by ``[referencename]_``, see, e.g., [Wall99]_.

- **Tables** are often entered in a rather strict format like this

  ::

     +----------+----------+
     + Header 1 | Header 2 |
     +==========+==========+
     | Left col | Right col|
     +----------+----------+

  However, this method is very time-consuming, and it may be difficult to add a new row or column,
  since the vertical lines must be perfectly aligned. The recommended since easier and more versatile way
  of adding a table is the following:

  ::

     .. list-table::
     :header-rows: 1

     * - Header 1
       - Header 2
     * - Left col
       - Right col


- **Images** can be entered in the following way:

  ::

     .. figure:: /_assets/figure.jpg
        :alt: This is the alternative name, if the figure cannot be shown
        :width: with should be given in pixels (e.g., 400px), but can also be given in %.
        :align: "top", "middle", "bottom", "left", "center", or "right"

        After one empty line, a figure caption is given (like this one).
        Beware that the indentation must not change!

  Note that figures cannot natively be entered in markdown. However, there is a way to enter them anyway by declaring a restructuredText element within the markdown file:

  .. code-block:: markdown

     ```{eval-rst}
     .. figure:: /_assets/figure.jpg
        :alt: This is the alternative name, if the figure cannot be shown
        :width: with should be given in pixels (e.g., 400px), but can also be given in %.
        :align: "top", "middle", "bottom", "left", "center", or "right"

        After one empty line, a figure caption is given (like this one).
        Beware that the indentation must not change!
     ```

  There may be other feature not included in markdown which can be entered in markdown by the exact same way.

- **Math** can be inserted either within the text or as separate equations.
  For inline math one may use the ``math`` rule like ``:math:`f(x) = x^2```.
  Separate equations are written as a directive ``.. math::``, for example:

  ::

     .. math::

        f(x) = \int_{\partial \Omega} \sigma_{ij}(\mathbf{x}) : \varepsilon_{ij}(\mathbf{x}) \mathrm{d} A

  which leads to

  .. math::

     f(x) = \int_{\partial \Omega} \sigma_{ij}(\mathbf{x}) : \varepsilon_{ij}(\mathbf{x}) \mathrm{d} A

- **Files for Download** should be included in the ``doc/documentation/src/assets`` directory.
  Within the text, you may use the sphinx ``:download:`` role, as shown :download:`here </_assets/testfile.pdf>`.
  This link has been created by ``:download:`here </_assets/testfile.pdf>```.

- **Need more html style?** If you need some special html coding, you may put this into ``<4C-sourcedir>/doc/documentation/src/_static/html_extensions.css``.
  For example, in the :ref:`coverage report <coveragereport>`, we needed green color, which comes in handy by adding

  ::

      .green {
          color: green;
      }

  to the css file.
  For using such a class in your text, you need to define a restructuredText role that uses this class,
  which you'll find in the definition of the string variable ``rst_prolog`` in ``<4C-sourcedir>/doc/documentation/conf.py.in``.
  Finally, you may use the newly defined role as ``:green:`this green phrase```, which gives you :green:`this green phrase`.
