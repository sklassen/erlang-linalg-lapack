-module(linalg_lapack_tests). 
-import(linalg_lapack,[transpose/1,matmul/2]).
-include_lib("eunit/include/eunit.hrl").

%transpose_1_test() ->
%	?assert(transpose([[8.0]])=:=[[8.0]]).
%
%transpose_2_test() ->
%	?assert(transpose([[1.0,2.0],[3.0,4.0]])==[[1.0,3.0],[2.0,4.0]]).
%
%transpose_3_test() ->
%	?assert(transpose([[1.0,2.0,3.0],[4.0,5.0,6.0]])==[[1.0,4.0],[2.0,5.0],[3.0,6.0]]).

multiply_2_test()->
	?assert(matmul([[1.0,2.0],[3.0,4.0]],[[1.0,3.0],[2.0,4.0]])==[[5.0,11.0],[11.0,25.0]]).

