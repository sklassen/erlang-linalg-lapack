-module(linalg_lapack).
-export([eye/1,dot/2,transpose/1,inv/1,matmul/2,svd/1]).

-on_load(init/0).


init() ->
    Directory=filename:dirname(code:which(linalg_lapack)),
	erlang:load_nif(Directory++"/../priv/linalg_lapack_nif", 0).

eye(_) -> 
        exit(nif_library_not_loaded).

dot(_,_) -> 
        exit(nif_library_not_loaded).

transpose(_) -> 
        exit(nif_library_not_loaded).

inv(_) -> 
        exit(nif_library_not_loaded).

matmul(_,_) -> 
        exit(nif_library_not_loaded).

svd(_) -> 
        exit(nif_library_not_loaded).


