-module(matrix_test).
-vsn('2011.1201').
-author('simon.klassen@adyne.com').

-export([start/0]).

%% s$u %*% diag(s$d) %*% t(s$v)


start() ->
	 X = [[7.52, -0.76,  5.13, -4.75,  1.33, -2.40] ,[-1.10,  0.62,  6.62,  8.52,  4.91, -6.77] ,[-7.95,  9.34, -5.66,  5.75, -5.49,  2.34] ,[1.08, -7.10,  0.87,  5.30, -3.52,  3.95]],

	 {{d,D},{u,U},{v,V}}= lapack:svd(X),


	 matrix:mmultiply(matrix:mmultiply(U,matrix:diag(D)),matrix:transpose(V)).
	 lapack:mmultiply(lapack:mmultiply(U,matrix:diag(D)),matrix:transpose(V)).

	c(matrix).
	c(matrix, [native]).
	c(matrix,[{hipe,o1}]).


	 lapack:mmultiply([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).
	 matrix:mmultiply([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).
	 matrix:mmultiply2([[1.0,2.0],[3.0,4.0]],[[1.0,2.0,3.0],[4.0,5.0,6.0]]).

	 matrix:mmultiply(matrix:sequential(2,2),matrix:sequential(2,3)).

	 matrix:mmultiply(matrix:sequential(200,200),matrix:sequential(200,400)).
	 lapack:mmultiply(matrix:sequential(200,200),matrix:sequential(200,400)).

	 element(1,timer:tc(matrix,mmultiply,[matrix:sequential(50,50),matrix:sequential(50,100)])).
	 element(1,timer:tc(lapack,mmultiply,[matrix:sequential(50,50),matrix:sequential(50,100)])).
	
	
