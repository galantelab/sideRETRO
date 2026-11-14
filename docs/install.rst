.. _chap_installation:

************
Installation
************

**sideRETRO** stores its source code on `github <https://github.com/galantelab/sideRETRO>`_
and uses `Meson build system <https://mesonbuild.com/>`_ to manage
configuration and compilation process.

Building requirements
=====================

- `Python 3 <https://www.python.org/>`_
- `Ninja <https://github.com/ninja-build/ninja/>`_

The building requirements for **Meson** can be obtained using package manager or
from source. For example, using `Ubuntu <https://ubuntu.com/>`_ distribution::

  $ sudo apt-get install python3 \
                         python3-pip \
                         python3-setuptools \
                         python3-wheel \
                         ninja-build

Installing Meson
================

**Meson** can be installed from the distribution packages::

  $ sudo apt install meson

or, for the most up-to-date version, through :code:`pip3`::

  $ pip3 install --user meson

In the latter case, remember to add :file:`~/.local/bin` to your
:code:`PATH`. For more information about using and installing Meson,
see: https://mesonbuild.com/Quick-guide.html

Project requirements
======================

- `zlib <https://www.zlib.net/>`_
- `HTSlib <http://www.htslib.org/>`_
- `SQLite3 <https://www.sqlite.org/>`_

If any requirements are not installed, during building, sideRETRO
will **download**, **compile** and statically link against the library.

Compiling and installing
========================

First, you need to **clone** sideRETRO repository::

  $ git clone https://github.com/galantelab/sideRETRO.git

Inside :file:`sideRETRO` folder, **configure** the project with Meson::

  $ meson setup build

if everything went well, you will see a new **folder** :file:`build`.
Now it is time to **compile** the code::

  $ meson compile -C build

It will be created the **executable** :file:`sider` inside the folder
:file:`build/src/`, which can already be used. Anyway, if want to
**install** sideRETRO to a system folder::

  $ sudo meson install -C build

By default, sideRETRO will install under :file:`/usr/local/`.
Choose another **prefix** when configuring the project::

  $ meson setup build --prefix=$HOME/.local

or change it, and see the other options, in an already configured
:file:`build` folder::

  $ meson configure build
  $ meson configure build --prefix=$HOME/.local

Testing
=======

The test suite has three parts: The **unit** tests, which need
`libcheck <https://libcheck.github.io/check/>`_, the **simulation**
tests, which run simulated retrocopies through the whole pipeline,
and the **cli** tests, which run each command with its options and
errors. The last two need only **Python 3**. Run them with::

  $ meson test -C build

or one suite at a time::

  $ meson test -C build --suite unit
  $ meson test -C build --suite simulation
  $ meson test -C build --suite cli

See :ref:`A quick test with simulated data <sim_test>` for running the
simulator by hand.
