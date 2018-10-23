-module(matrix_lapack).
-export([transpose/1,multiply/2,svd/1]).

-on_load(init/0).


init() ->
    Directory=filename:dirname(code:which(matrix_lapack)),
	erlang:load_nif(Directory++"/../priv/matrix_lapack_nif", 0).

transpose(_) -> exit(nif_library_not_loaded).
multiply(_,_) -> exit(nif_library_not_loaded).
svd(_) -> exit(nif_library_not_loaded).


