# Roadmap for the trusat-orbit repository code:

## Near Term (next few weeks)
1. [x] Pull out loop-intensive functions into accelerated module
1. [ ] Set up Travis Continuous Integration Testing
1. [ ] Publish pypi trusat (trusat-orbit?) module
1. [ ] Expand documentation for using trusat.tle and trusat.iod modules
1. [ ] Create jupyter notebook for development, experimentation, visualization
1. [ ] Restore functionality of local IOD file processing with satfit. [Issue #8](https://github.com/consensys-space/trusat-orbit/issues/8)
1. [ ] Fix satid/elfind compatibility with python-skyfield
1. [ ] Incorporate XF functionality into TLE processing scripts (See tle_util.py header comments)
1. [ ] Transition python-skyfield dependencies to native python-sgp4 calls


## Mid Term (next few months)
1. [ ] Extend satfit functionality in PyQT module
1. [ ] Write verison of TLE Retrieve in PyQT, utilizing trusat.tle_util
1. [ ] Merge trusat-orbit and trusat-backend into a common repo and package
1. [ ] API access to all data (after SeeSat-L user opt-in/opt-out)
1. [ ] Secure compute for SITE data - allow geocentric (ICRS) observation vectors / observer location


## Long Term
1. [ ] Transition to modules compatible with OREKit and [UT-Austin orbdetpy](https://github.com/ut-astria/orbdetpy)
1. [ ] Full automation with user rank/trust, object priority and object confidence
1. [ ] Parallelized processing for super-compute applications
