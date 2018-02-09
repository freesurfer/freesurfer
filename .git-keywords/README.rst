Helper scripts for keyword expansion in git
###########################################

:date: 2015-05-14
:tags: git, keywords
:author: Roland Smith

.. Last modified: 2015-05-14 18:02:25 +0200

One of the things I liked about the old rcs_ revision control system was that
it supported keyword expansion in files.  Unlike systems like ``rcs``, ``cvs``
and ``subversion``, the ``git`` revision control system cannot provide keyword
expansion. The cause for this is that you can't modify a file with information
about the commit after you've committed, because ``git`` checksums the file
first.

.. _rcs: http://en.wikipedia.org/wiki/Revision_Control_System

Git will let you inject text in a file when it is checked out, and remove it
when it is checked in. There are two ways of doing this. First, you can use
the ``ident`` attribute_. For any file type that has the ``ident`` attribute
set (in ``.gitattributes``), git will look for the string ``$Id$`` on checkout
and add the SHA-1 of the blob to it like this: ``$Id:
daf7affdeadc31cbcf8689f2ac5fcb6ecb6fd85e $``. While this unambiguously
identifies the commit, it is not all that practical.

* It cannot tell you the relative order of two commits.
* It doesn't tell you the commit date.

Luckily, keyword expansion can be done with git using attributes_.

.. _attribute: http://git-scm.com/book/en/v2/Customizing-Git-Git-Attributes
.. _attributes: http://git-scm.com/book/en/v2/Customizing-Git-Git-Attributes

In my global git configuration file (``~/.gitconfig``) I have defined a
filter_ called "kw":

.. _filter: http://git-scm.com/docs/gitattributes

.. code-block:: ini

  [filter "kw"]
     clean = kwclean
     smudge = kwset

This configuration uses two programs (which should be in your ``\$PATH``)
called ``kwset`` and ``kwclean`` to expand and contract keywords. These are
two scripts written in python_ 3.

.. _python: http://python.org/

To *enable* these substitutions, you have to use git attributes. E.g. to have
keyword substitutions in *all* files in a repository, you need to add the
following to the ``.gitattributes`` file in that repository;

.. code-block:: ini

    * filter=kw

Such a general use of filters can be problematic with e.g. binary files like
pictures. As a rule, modifying the contents of a binary (especially adding or
removing bytes) tends to *break* them.

It is therefore better to be explicit and specific as to what types of file
the filter should apply to;

.. code-block:: ini

    *.py filter=kw
    *.txt filter=kw

With this filter setup, file types that contain keywords and which are listed
as such in the ``.gitattributes`` file will have them expanded on checkout.

To make these updated keywords visible in the working directory, changed
objects will have to be checked out after their changes have been committed.
To accomplish this, we can use the ``post-commit`` hook_. There are several
possible choices here. You can e.g.:

.. _hook: http://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks

* Check out the files which have changed since the previous commit.
* Check out *all* files.

The first one is probably the most common case. I wrote the script
``update-modified-keywords.py`` for it. After a check-in it checks out all the
files that were modified in the last commit.

But if all the directories in one file are part of one project, you probably
want all files to carry the same date/revision. This is what the
``update-all-keywords.py`` script is for. After a check-in it checks out all the
files that are under git's control.

Put both these scripts in a location in your ``$PATH``, and then make symbolic
links from ``.git/hooks/post-commit`` to the appropriate script.

.. NOTE::

  .. image:: http://i.creativecommons.org/p/zero/1.0/88x31.png
        :alt: CC0
        :align: center
        :target: http://creativecommons.org/publicdomain/zero/1.0/

  To the extent possible under law, Roland Smith has waived all copyright and
  related or neighboring rights to ``kwset.py``, ``kwclean.py``,
  ``update-all-keywords.py`` and ``update-modified-keywords.py``. These
  works are published from the Netherlands.
