# Finite Element libraries

## Description
This project contend example of program that illustrated use of :

-  FEniCSx (0.8.0)
-  MFEM (4.7.0)

A document (in tex format) collect comparison of those libraries. 

## Installation
setting.mk.in must be copied into setting.mk in each directory and adapted to your own ecosystem  (path,dependency,...).
There are many macros at the beginning of the source files that need to be set to tune the tests.
Then make launched in root folder should normally compile examples. 

## Usage
Executable can be launched in each directory containing a test case.


## Support
to be discuss

## Roadmap
Add examples to test topics not already illustrated by current version.

## Contributing
feel free

## License
Copyright (c) 2025      - Ecole Centrale de Nantes

Doc (latex file) is under CC-BY-NC-ND <img src="doc.data/img/cc-by-nc-nd.png" alt="CC-BY-NC-ND" width="50" height="18">



All c++/python sources are under GNU LGPL license given in LICENCE.md

Except admfem.hpp/taddensemat.hpp/tadvector.hpp which are under BSD 3-Clause License with the following Copyright:
Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC All rights reserved.
They have been copied to insure that Automatic Differentiation test with MFEM can be done without searching for those files 
not provided yet in the core source file folder of the libraries.

## Project status
early birth:
It is the result of internal work by some members of the MECNUM research team at the GeM Institute (https://gem.ec-nantes.fr), who wish to provide researchers and PhD students in the team with some guidance on these libraries. It is made available to the public through this website for use both internally and potentially externally.

