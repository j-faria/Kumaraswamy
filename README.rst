Kumaraswamy
===========

A simple implementation of the Kumaraswamy distribution.

|License MIT| |PyPI version| |DOIzenodo|



How to use
----------

Install it from pip (**kumaraswamy** only depends on numpy)

::

    pip install kumaraswamy

and it's ready to use from Python

.. code:: python

    import kumaraswamy

The package provides one simple class called ``kumaraswamy``, which implements the distribution.
It is intended to mimic the API of ``scipy.stats``.

.. code:: python

    from kumaraswamy import kumaraswamy

    d1 = kumaraswamy(a=0.5, b=0.5)
    
the ``d1`` object now has methods

-  ``pdf(x)`` and ``logpdf(x)``
-  ``cdf(x)``
-  ``rvs(size)``

to calculate the probability density function (and its logarithm), the
cumulative density function, and to get random samples from the
distribution.

Also available are some basic properties specific to the distribution

.. code:: python
    
    d1.mean
    d1.var
    d2.mode  # see help(d2.mode) for details



License
-------

Copyright 2018 João Faria.

**kumaraswamy** is free software made available under the MIT License. For
details see the LICENSE_ file.

.. _License: https://github.com/j-faria/Kumaraswamy/blob/master/LICENSE
.. |License MIT| image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/j-faria/Kumaraswamy/blob/master/LICENSE
.. |PyPI version| image:: https://badge.fury.io/py/kumaraswamy.svg
   :target: https://pypi.org/project/Kumaraswamy/
.. |DOIzenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3951068.svg
   :target: https://doi.org/10.5281/zenodo.3951068
