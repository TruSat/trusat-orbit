![TruSat banner](https://trusat-assets.s3.amazonaws.com/readme-banner.jpg)

# trusat-orbit

This repo supports the analyses of [IOD/RDE/UK positional formatting formats](http://www.satobs.org/position/posn_formats.html) and generation of TLEs at https://TruSat.org/

Currently, its orbit propagation code is based on a Python port of [Scott Campbell's C++ satfit code base]( https://github.com/interplanetarychris/scottcampbell-satfit). After initial prototyping, it is an aim of this repo to include OREKit and related tools for more advanced processing of orbit-related calculations.

# Getting started with TruSat-orbit
First, we recommend setting up a [python virtual environment](https://realpython.com/python-virtual-environments-a-primer/)
```
git clone https://github.com/TruSat/trusat-orbit.git
cd trusat-orbit
pip install -r requirements.txt
python trusat/satfit.py
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
* tests_profile.py

Run with python -m unittest discover tests/

## Dependencies
* Uses Brandon Rhodes [python-sgp4](https://github.com/brandon-rhodes/python-sgp4) with C++ accelerations
* (Currently) requires Cython for the c-accelerated analyses module
* (Currently) assumes a connection to a database, see [this PR] (https://github.com/TruSat/trusat-backend/pull/111) for setting up your own local copy
* (Currently) requires pulling [trusat-backend](https://github.com/TruSat/trusat-backend) into the same parent directory as trusat-orbit

## Coding Style
Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for any Python code and the style guide recommended for any other language.

Additionally see [Best of the Best Practices](https://gist.github.com/sloria/7001839)

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
