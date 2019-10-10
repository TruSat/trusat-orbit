# trusat-orbit

This repo supports the analyses of [IOD/RDE/UK positional formatting formats](http://www.satobs.org/position/posn_formats.html) and generation of TLEs at http://TruSat.org/

Currently, its orbit determination code is based on a Python port of [Scott Campbell's C++ satfit code base]( https://github.com/interplanetarychris/scottcampbell-satfit). After initial prototyping, it is an aim of this repo to include OREKit and related tools for more advanced processing of orbit-related calculations.

# Getting started with TruSat-orbit
First, we recommend setting up a [python virtual environment](https://realpython.com/python-virtual-environments-a-primer/)
```
git clone https://github.com/consensys-space/trusat-orbit.git 
pip install -r requirements.txt
python satfit.py
```

## Contents
* **iod.py** - Utilities for importing, validating, and operating on IOD/RDE/UK positional formatting formats 
* **tle_util.py** - Utilities to import, export, validate and operate on Two-Line Element sets
* **satfit.py** - Suite of utilities based on and extending [Scott Campbell's C++ satfit code base]( https://github.com/interplanetarychris/scottcampbell-satfit) for reading visual observations and updating TLEs
  * **satid.py** - Search TLE catalog for possible match to an UNIDentified satellite TLE
  * **elfind.py** - Generate a provisional TLE from 2-3 IOD records

### Tests - Unit tests for the above
* tests_iod.py 
* tests_satfit.py
* tests_TLE.py

## Dependencies
This codebase incorporates an updated and cython-accelerated version of python-SGP4, from the branch [python-sgp4/cython-7-dec-15-vallado](https://github.com/interplanetarychris/python-sgp4/tree/cython-7-dec-15-vallado)

## Coding Style
Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for any Python code and the style guide recommended for any other language.

## Maintaining Repo
[Style Guide](https://github.com/agis/git-style-guide)
With the addition of commits to the master branch are done through PRs (Pull Request).
## Releasing Versions
1. Checkout master
2. pull from repo
3. run the unittests
4. create a tag with the new version number, starting with a 'v'. eg:

```git tag v0.1.1 -m "Version 0.1.1```
[Version Numbering](semver.org)
5. push changes to github `git push --follow-tags`
7. check verification tools

## Running tests

There are two separate ways to use `test_TLE.py`. They do very different things.

To run the tests, use Python3's `unittest` module. E.g.: 

`python3.7 -m unittest -v test_TLE.Tests`

To execute the `main()` function, just execute it. E.g.:

`python3.7 test_TLE.py`

!TODO: what's the reason to do the latter? Why doesn't the latter do the same as the former?

### Include database tests

Set the following environment variables before running your test:

- `TRUSAT_DATABASE_NAME`
- `TRUSAT_DATABASE_HOST`
- `TRUSAT_DATABASE_USER`
- `TRUSAT_DATABASE_PASSWORD`

E.g. run:

`TRUSAT_DATABASE_NAME=opensatcat_local TRUSAT_DATABASE_HOST=127.0.0.1 TRUSAT_DATABASE_USER=test TRUSAT_DATABASE_PASSWORD=test python3.7 -m unittest -v test_TLE.Tests`