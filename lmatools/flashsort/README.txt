The "autosort" directory contains support code for flash sorting conducted by a FORTRAN algorithm, which is not included due to license ambiguity for that algorithm.

The only missing file is mflash.py, which contains a long string with the fortran flash-sorting algorithm. Many parameters are encoded as %(param_name)s so that a string format operation can replace those parameters and write out mflash.f, which is then compiled and run by autorun.py. As currently configured, the FORTRAN routine takes input from stdin, which is marshaled by code in autorun.py

This code is licensed under a BSD-type license.

-ECB, 21 Nov 2010