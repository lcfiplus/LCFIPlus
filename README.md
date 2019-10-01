# LCFIPlus
[![Build Status](https://travis-ci.org/lcfiplus/LCFIPlus.svg?branch=master)](https://travis-ci.org/lcfiplus/LCFIPlus)
[![Build Status](https://scan.coverity.com/projects/14336/badge.svg)](https://scan.coverity.com/projects/lcfiplus-lcfiplus)

Flavor tagging code for ILC detectors, for documentation consult confluence at [https://confluence.slac.stanford.edu/display/ilc/LCFIPlus](https://confluence.slac.stanford.edu/display/ilc/LCFIPlus)

LCFIPlus is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## License and Copyright
Copyright (C) 2005-2017, LCFIPlus Authors

LCFIPlus is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.


## Release notes

see: [./doc/ReleaseNotes.md](./doc/ReleaseNotes.md)

Primary Vertices are built out of Tracks, then secondary vertices are found, in a further step vertices and jets are matched to each other to find classification probabilities that the jet originates from a b-quark, or a c-quark, or a light flavor quark

JetFinderFix Version of the code evaluates methods which mix slightly different input collection for vertexing and jet clustering. The original code expected the inputs to be the same, thus it had to be modified to avoid double counting of inputs.
