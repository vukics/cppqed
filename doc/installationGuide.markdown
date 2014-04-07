Installation guide {#installationguide}
==================

# Quick start: Installation from packages (Ubuntu, Arch Linux) {#installationguidequickstart}

Binary packages have been prepared for Arch Linux and [Ubuntu](http://www.launchpad.net/~raimar-sandner/cppqed),
which might also work on Debian (not tested).  Installation from binary packages is recommended for users
who want to try C++QED or use it exclusively on its highest level, i.e. writing
scripts which use elements and interactions already implemented in C++QED. If
you are not using Ubuntu or want to develop your own elements and interactions
you have to install C++QED from source.

### Ubuntu

All packages are compiled for Ubuntu Precise (12.04 LTS) and Saucy (13.10). You can choose between the latest release and daily builds. For the release:

    $ sudo add-apt-repository ppa:raimar-sandner/cppqed-development
    $ sudo apt-get update
    $ sudo apt-get install libcppqedelements-2.99-dev cppqedscripts python-cpypyqed

For the daily builds:

    $ sudo add-apt-repository ppa:raimar-sandner/cppqed-daily
    $ sudo apt-get update
    $ sudo apt-get install libcppqedelements-daily-dev cppqedscripts python-cpypyqed

The current documentation can be installed with

    $ sudo apt-get install cppqed-doc

### Arch

The C++QED package is in the [Arch User Repository (AUR)](https://aur.archlinux.org/packages/cppqed-git), you can install it for example with the `yaourt` package manager (`pacman` alternative which can install from AUR):

    $ yaourt -S cppqed-git

This will install all dependencies and compile C++QED.

# Installation from source {#installationguidefromsource}

Foo Bar