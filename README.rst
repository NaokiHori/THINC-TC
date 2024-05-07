########
THINC-TC
########

|License| |LastCommit| |CI| |DOCS|

.. |License| image:: https://img.shields.io/github/license/NaokiHori/THINC-TC
.. _License: https://opensource.org/licenses/MIT

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/THINC-TC/main
.. _LastCommit: https://github.com/NaokiHori/THINC-TC/commits/main

.. |CI| image:: https://github.com/NaokiHori/THINC-TC/actions/workflows/ci.yml/badge.svg

.. |DOCS| image:: https://github.com/NaokiHori/THINC-TC/actions/workflows/documentation.yml/badge.svg
.. _DOCS: https://naokihori.github.io/THINC-TC

.. image:: https://github.com/NaokiHori/THINC-TC/blob/main/docs/source/thumbnail.jpg
   :target: https://youtube.com/shorts/GlfF6iDPeiw
   :width: 100%

********
Overview
********

A simplified / polished implementation of the THINC (volume-of-fluid) method for Taylor-Couette flows, which is briefly discussed in the appendix of `this manuscript <https://doi.org/10.1017/jfm.2023.29>`_.

Basically this is the combination of `the single-phase Taylor-Couette solver <https://github.com/NaokiHori/SimpleTCSolver>`_ and `the volume-of-fluid solver for Cartesian domains <https://github.com/NaokiHori/SimpleVOFSolver>`_.
Its aim is to perform interface-resolved direct numerical simulations of two-liquid turbulent Taylor-Couette flows without density contrast.

***********
Quick start
***********

Check the documentation of the aforementioned projects.

