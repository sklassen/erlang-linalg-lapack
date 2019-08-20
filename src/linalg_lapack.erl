-module(linalg_lapack).
-export([transpose/1,matmul/2,svd/1]).

-on_load(init/0).


init() ->
    Directory=filename:dirname(code:which(linalg_lapack)),
	erlang:load_nif(Directory++"/../priv/linalg_lapack_nif", 0).

transpose(_) -> exit(nif_library_not_loaded).
matmul(_,_) -> exit(nif_library_not_loaded).
svd(_) -> exit(nif_library_not_loaded).


