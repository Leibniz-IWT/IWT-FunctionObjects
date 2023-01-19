# IWT-FunctionObjects
IWT OpenFOAM FunctionObjects Library

## Contributing Instructions
1. Create a new folder with a clear name describing what your functionObject does.
2. Edit the Make/files file adding your .C file to be compiled in the list.
3. Edit the Make/options file adding any libraries necessary for your new functionObject. Make sure all the libraries you add are absolutely necessary!
4. Edit Allwmake adding symbolic links for your .C and .H files to lnInclude.
5. Run Allwmake to check your and all preexisting functionObjects compile without errors.
6. Run Allwclean to remove the compiled binaries and the lnInclude folder.
7. Submit your changes in a new branch in the repository and create a pull request to be included in the next main branch's next release.

## Installation Instructions
### Prerequisites
- [OpenFOAM v6](https://openfoam.org/download/source/ "OpenFOAM v6") installed
- [OpenFOAM v6](https://openfoam.org/download/source/ "OpenFOAM v6") sourced in the currently open shell

### Downloading
- `git clone https://github.com/Leibniz-IWT/IWT-FunctionObjects.git`

or
- Manually download and extract the release zip archive

### Compilling
1. `cd IWT-FunctionObjects`
2. `./Allwmake`
