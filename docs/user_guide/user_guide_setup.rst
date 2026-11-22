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
The Ring-Moon Systems Node publishes these resources, the kernels and ``SPICE.db``
included, at ``https://storage.googleapis.com/rms-node-oops-resources``, a URL that
``OOPS_RESOURCES`` can name directly.

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
defines the camera's fields of view. To choose the bodies that OOPS defines, to use
the kernels released before a given date, or to furnish a metakernel instead of the kernels
chosen from the database, call
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.initialize` before opening the
first image. Its ``planets`` option takes one planet number, from 1 to 9, or a sequence
of them. OOPS always defines the Sun, Mercury, Venus, Earth, and the Moon, and this call
adds Jupiter and Pluto, planets 5 and 9, with their moons:

.. code-block:: python

   from host_new_horizons import NewHorizonsLORRI

   NewHorizonsLORRI.initialize(planets=(5, 9))

Its ``asof`` option takes a date, such as ``'2016-01-01'``, and each kernel is then
chosen from those released before that date. Where every kernel that a query could
choose was released on or after the date, the query silently uses the kernel it would
choose without a date. That happens only when the whole query comes back empty: a set
of SP or C kernels that covers the mission keeps just the kernels released before the
date, so a date inside the mission can leave the later part of the mission uncovered.
With ``asof='2015-06-01'``, for example, the C kernels end on 2015-01-01. The C kernels
never fall back to the newest: for any date up to and including 2015-03-12, the day
the first of them was released,
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.initialize` raises
:exc:`ValueError`. The date decides the leap seconds too: with ``asof='2015-06-01'``,
the leap seconds kernel released on 2015-01-05 is used, which lacks the leap second at
the end of 2016.

With ``meta``, the date chooses the metakernel when ``meta`` names one in the kernel
database: the newest one released before the date is used, and if none was,
:meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.initialize` raises
:exc:`ValueError`, because this choice never falls back to the newest either. A
metakernel given by its file path is used whatever the date. The kernels that a
metakernel lists are furnished whatever their release dates, and a leap seconds kernel
among them takes the place of the one that the date chose.

Later calls to :meth:`~host_new_horizons.oops.lorri.NewHorizonsLORRI.initialize` leave
the first call's definitions and kernels as they are.

Image files and indexes
=======================

The archive itself is found without an environment variable. Each FITS file, PDS3
label, or index is passed by its path or URL, and a remote file is retrieved and cached
when it is read.
