The Knockoff Filter for MATLAB
==========================

This package provides a versatile MATLAB interface to the knockoff methodology.

# Installation

## Requirements

Basic software requirements:

- MATLAB 2012 or later
- Statistical toolbox
- CVX (optional; required for creating SDP knockoffs)

Warning: your CVX release must be from April 2014 or later. Earlier
versions of CVX contain a serious bug.
  
Additional software requirements:

 - glmnet for MATLAB
 
Many statistics used by the knockoff filter are computed with glmnet, 
which can be downloaded from:
http://web.stanford.edu/~hastie/glmnet_matlab/download.html

If you have not already installed the glmnet package for MATLAB, please download it using the link provided above 
and unpack the archive into the current directory. 
This should create a subdirectory named "glmnet_matlab".

## Installation from source

You can install the lastest development version by cloning this repository, or downloading it as an archive, and then loading the source code from MATLAB.

To install the package, save the directory "knockoffs_matlab" of this repository anywhere on your machine, say to `<path>`, then add the lines
```Matlab
addpath("<path>/knockoffs_matlab")
addpath("<path>/knockoffs_matlab/glmnet_matlab")
```

To test your installation, you can run one of the demo files in the "examples" subdirectory. Alternatively, if you have MATLAB 2013a or newer, you can execute the test suite by typing
```Matlab
runtests('knockoffs.tests')
```

in your MATLAB command window. All the tests should pass.

## Documentation

For an overview of the functions in this package, run

```Matlab
help knockoffs
```

For a list of included knockoff statistics, run

```Matlab
help knockoffs.stats
```

Besides the documentation associated with individual functions 
(accessible via the "help" function), the main source of documentation
is the collection of examples in the "examples" subdirectory.

The file "examples_basic.m" is a good place to start.

## Resources
For more information and tutorials, visit
https://web.stanford.edu/group/candes/knockoffs/software/knockoffs/download-matlab.html

## Credits

This package was developed by Matteo Sesia, Lucas Janson, Emmanuel Candès, Yingying Fan and Jinchi Lv.

An earlier version of this package was developed by Evan Patterson, based on code originally written Rina Foygel Barber, Emmanuel Candès and Evan Patterson: https://bitbucket.org/epatters/knockoff-filter.

## License

This software is distributed under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) and it comes with ABSOLUTELY NO WARRANTY.