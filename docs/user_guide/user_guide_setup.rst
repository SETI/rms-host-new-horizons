=====
Setup
=====

Installation
============

The package needs Python 3.11 or later and installs from PyPI:

.. code-block:: sh

   pip install rms-host-new-horizons

It brings in OOPS (``rms-oops``) and the readers that the host uses for PDS3 labels,
FITS files, and PDS3 index tables.

SPICE kernels
=============

The geometry of every image comes from the New Horizons SPICE kernels. The host reads
them through ``spicedb``, which is distributed with OOPS and finds them with two
environment variables:

``SPICE_PATH``
    The directory that holds the kernels, as a local path or a URL. If it is not set,
    ``$OOPS_RESOURCES/SPICE`` is used.

``SPICE_SQLITE_DB_NAME``
    The SQLite database that lists the kernels. If it is not set, ``SPICE.db`` in the
    ``SPICE_PATH`` directory is used.

The kernels and the database are the ones that OOPS itself uses. The
`rms-oops README <https://github.com/SETI/rms-oops#environment-variables>`_ describes
these variables and the ``SPICE`` directory of ``$OOPS_RESOURCES`` that holds them.

The host selects its kernels from the database once, for the whole mission from the
launch on 2006 January 19 through 2021: the leap-second kernel, the planetary constants
of the Jupiter and Pluto systems, the spacecraft's instrument, frame, and clock kernels,
its predicted and reconstructed ephemerides, and its C kernels. A metakernel can be
furnished in their place (see below).

Registering the host
====================

The host classes import lazily. ``import host_new_horizons`` gives access to
:class:`~host_new_horizons.oops.lorri.NewHorizonsLORRI` without importing OOPS until the
class is first used. Importing ``host_new_horizons.oops`` registers the host with
:class:`~oops.host.Host`, so that :meth:`oops.host.Host.from_file` recognizes a LORRI
image:

.. code-block:: python

   import oops
   import host_new_horizons.oops

   obs = oops.Host.from_file(
       'NHJULO_2001/data/20070110_003071/lor_0030713591_0x633_sci.fit')

Initialization
==============

The first call to :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_file` or
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.from_index` loads the kernels and
defines the camera's fields of view. To choose the bodies that OOPS defines, to use only
the kernels that existed before a given date, or to furnish a metakernel instead of the
kernels chosen from the database, call
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.initialize` before opening the
first image. Its ``planets`` option takes one planet number, from 1 to 9, or a sequence
of them. This call defines Jupiter and Pluto, planets 5 and 9, and their moons, and no
other planet:

.. code-block:: python

   from host_new_horizons import NewHorizonsLORRI

   NewHorizonsLORRI.initialize(planets=(5, 9))

Later calls to :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.initialize` have no
effect.

Image files and indexes
=======================

The host needs no environment variable for the archive itself. Each FITS file, PDS3
label, or index is passed by its path or URL, and a remote file is retrieved and cached
when it is read.
