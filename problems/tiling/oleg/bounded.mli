(* Partially known dynamic values for safe array indexing *)


type bound			(* A single bound, abstract *)

val string_of_bound : bound -> string

				(* Partial order on bounds *)
type fuzzy_bool = True | False | Undecided
val leq_bound : bound -> bound -> fuzzy_bool

val offset_bound : int -> bound -> bound

type bounds = bound * bound	(* Lower and higher bounds *)

val string_of_bounds : bounds -> string

val array_bounds : bounds	(* Full bounds of an array *)

				(* A bounded dynamic value of the type int *)
				(* Hereafter, 'c stands for the environment 
				   classifier. *)
type 'c bounded

val string_of_bounded : 'c bounded -> string

				(* Extract the dynamic value from bounded *)
val dyn_of_bounded : 'c bounded -> ('c,int) code

				(* Compare bounded: we compare bounds 
				   rather than code values *)
val same_bounded : 'c bounded -> 'c bounded -> bool

				(* Check to see if one bounded immediately
				   precedes the other. That is, the upper 
				   bound of one matches with the lower bound
				   of the other. *)
val successive_bounded : 'c bounded -> 'c bounded -> bool

val offset_bounded : int -> 'c bounded -> 'c bounded

				(* The bounded value corresponding to the
				   lower bound of an array *)
val min_bounded : 'c bounded

(* Check to see if a given bounded value is within the given bounds. 
   Return:
     - None if it is definitely out of bounds
     - Some c where c is the code (dynamic value) if it is definitely
       within the bounds
     - raise an error if cannot decide
*)
val within_bounds : 'c bounded -> bounds -> ('c,int) code option

				(* Wrapper that produces the top-level 
				   function for the stencil computation with 
				   one input and one output array,
				   given the generator for the computation *)
val make_fun : 
  ('c bounded -> 			(* for the upper bound of arrays *)
   ('c,float array) code -> 		(* input array *)
   ('c,float array) code ->		(* output array *)
   ('c,'w) code) ->
  ('c, int -> float array -> float array -> 'w) code

				(* For-loop with bounded index bounds *)
val forloop : 'c bounded -> 'c bounded -> ('c bounded -> ('c,unit) code)
  -> ('c,unit) code
