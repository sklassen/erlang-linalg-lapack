lapack - NIFs for Erlang
============================

This is an erlang matrix library that uses the LAPACK — Linear Algebra PACKage. LAPACK in turn uses blas and atlas. 

Installation
-----

This code was compiled and run under ubuntu. You will need blas and atlas development libraries.

	> sudo apt-get install libatlas-base-dev libblas-dev

Then, in the top directory, compile using rebar

	> rebar compile

Usage
-----

Only two functions are supported as yet.

	Erlang R16B03 (erts-5.10.4) [source] [64-bit] 

	Eshell V5.10.4  (abort with ^G)
	1> lapack:mmultiply([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).

	2>  X = [[7.52, -0.76,  5.13, -4.75,  1.33, -2.40] ,[-1.10,  0.62,  6.62,  8.52,  4.91, -6.77] ,[-7.95,  9.34, -5.66,  5.75, -5.49,  2.34] ,[1.08, -7.10,  0.87,  5.30, -3.52,  3.95]],
	3>  {{d,D},{u,U},{v,V}}= lapack:svd(X).


