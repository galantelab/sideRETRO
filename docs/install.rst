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

The recommended way to install the most up-to-date version of
**Meson** is through :code:`pip3`::

  $ sudo apt install meson
  
Or "$ pip3 install --user meson" But in this case, remember to set the environment variables. 

For more information about using and installing Meson, see:
https://mesonbuild.com/Quick-guide.html

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

  $ meson build

if everything went well, you will see a new **folder** :file:`build`.
Now it is time to **compiling** the code::

  $ ninja -C build

It will be created the **executable** :file:`sider` inside the folder
:file:`build/src/`, which can already be used. Anyway, if want to
**install** sideRETRO to a system folder::

  $ sudo ninja -C build install

By default, sideRETRO will install under :file:`/usr/local/`.
The configure script can **customize** the prefix directory. Run::

  $ meson build configure

for instructions and other installation options.
