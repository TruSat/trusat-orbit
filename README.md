# sathunt-iod
## Coding Style
Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for any Python code and the style guide recommended for any other language.
## Maintaining Repo
[Style Guide](https://github.com/agis/git-style-guide)
With the addition of commits to the master branch are done through PRs (Pull Request).
## Releasing Versions
Modified from [pyorbital](https://github.com/pytroll/pyorbital/blob/master/RELEASING.md)
1. checkout master
2. pull from repo
3. run the unittests
4. create a tag with the new version number, starting with a 'v'. eg:

```git tag v0.1.1 -m "Version 0.1.1```
[Version Numbering](semver.org)
5. push changes to github `git push --follow-tags`
7. check verification tools
