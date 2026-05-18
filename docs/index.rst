.. KN1D documentation master file, created by
   sphinx-quickstart on Fri Apr 17 10:35:42 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

KN1D documentation
==================

KN1D is a 1D-space, 2D-velocity kinetic neutrals code developed by B. LaBombard
(MIT PSFC). KN1DPy is a Python translation of the original IDL code (the
original IDL version is also included in the ``IDL/`` directory of this
repository).

This translation was written by Nicholas Brown at William & Mary
(njbrown@wm.edu), and is now maintained by Jamie Dunsmore at MIT
(jduns@mit.edu).

A comprehensive description of the algorithm can be found in this `PSFC report
(LaBombard, 2001) <https://library.psfc.mit.edu/catalog/reports/2000/01rr/01rr003/01rr003_full.pdf>`_.


.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   kn1d_lite

.. toctree::
   :maxdepth: 2
   :caption: Reference

   API Reference <api>
