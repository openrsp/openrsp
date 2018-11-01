.. _chapter_sphinx:

How Sphinx works
================

These pages are generated using Sphinx. If you want to find out more about
RST/Sphinx, please read http://sphinx-doc.org/rest.html. RST is a subset of
Sphinx. Sphinx is RST with some extensions.


How to modify the webpages
--------------------------

You can modify these pages by cloning our public documentation repository::

  git clone git@gitlab.com:openrsp/website.git

Once you commit and push, a post-receive hook
updates the documentation on http://openrsp.readthedocs.org.
This typically takes less than a minute.
Our main page http://openrsp.org redirects to http://openrsp.readthedocs.org.

How to locally test changes
---------------------------

You don't have to push to see and test your changes.
You can test them locally.
For this install python-sphinx and python-matplotlib.
Then build the pages with::

  make html

Then point your browser to _build/html/index.html.
The style is not the same but the content is what you
would see after the git push.
