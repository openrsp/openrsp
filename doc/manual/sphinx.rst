

How Sphinx works
================

These pages are generated using Sphinx. If you want to find out more about
RST/Sphinx, please read http://sphinx-doc.org/rest.html. RST is a subset of
Sphinx. Sphinx is RST with some extensions.


How to modify the website
-------------------------

The website is generated from RST sources under ``doc/``.
Once a pull request is merged, a post-receive hook
updates the documentation on https://openrsp.readthedocs.io.
This typically takes less than a minute.
Our main page http://openrsp.org redirects to https://openrsp.readthedocs.io.


How to locally test changes
---------------------------

You don't have to push to see and test your changes.
You can test them locally.
For this install the Python packages ``sphinx`` and ``sphinx_rtd_theme``.
Then build the pages with::

  $ sphinx-build doc/ _build

Then point your browser to ``_build/html/index.html``.
The style is not the same but the content is what you
would see after a successful pull request merge.
