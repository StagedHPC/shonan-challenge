(* Partially known dynamic values for array indexing

To make sure that the generated code has no array bound errors, we
need at the code generation time some information about dynamic values
that are to be used as indices of an array.  Specifically, we need the
lower (lwb) and the upper (upb) bounds of the value.  We always assume
that lwb < upb and the value v with the bounds (lwb,upb) satisfies
     lwb <= b < upb

In our application -- stencil computations -- all arrays have the same
length and are 0-indexed (that is, the smallest valid index into the
array is zero). The length of the array is not statically known. We do
assume that it is at least min_high (which we must check, see make_fun
below).

Bounds are treated as abstract outside this module. The constant
min_bounded gives the value corresponding to the smallest index
of the array, and wrappers like make_fun produce the corresponding
bounded value for the upper bound of the array. The loop generators
like forloop annotate the loop variable with the loop bounds.
*)

(* Configuration data *)

let min_high = 4;;	(* Minimal upper bound for an array    *)

(*
   The bounds themselves come in the form of an offset from
   min array bound and max array bound.
*)
type bound = 
  | Low of int  		(* offset *)
  | High of int
;;

type fuzzy_bool = True | False | Undecided;;

(* Partial order on bounds *)
let leq_bound : bound -> bound -> fuzzy_bool = fun lb ub ->
 match (lb,ub) with
 | (Low x, Low y)   -> if x <= y then True else False
 | (High x, High y) -> if x <= y then True else False
 | (Low x, High y)  -> if x <= min_high + y then True else Undecided
 | (High x, Low y)  -> if min_high + x > y then False else Undecided
;;

let offset_bound : int -> bound -> bound = fun n -> function
  | Low x  -> Low (x+n)
  | High x -> High (x+n)
;;

let string_of_bound : bound -> string = function
  | Low 0 -> "L"
  | Low n when n > 0 -> "L+" ^ string_of_int n
  | Low n when n < 0 -> "L-" ^ string_of_int (-n)
  | High 0 -> "H"
  | High n when n > 0 -> "H+" ^ string_of_int n
  | High n when n < 0 -> "H-" ^ string_of_int (-n)
;;


type bounds = bound * bound;;		(* Lower and higher bounds *)

let string_of_bounds : bounds -> string = fun (lwb,upb) ->
  "[" ^ string_of_bound lwb ^ "," ^ string_of_bound upb ^ ")"
;;


(* A bounded dynamic value of the type int *)
(* Hereafter, the type variable 'c stands for the environment classifier. *)
type 'c bounded = Bounded of
                    bounds *		(* lower and the upper bound *)
                    ('c,int) code 	(* the value itself *)
;;


let offset_bounded : int -> 'c bounded -> 'c bounded = function
  | 0 -> fun bv -> bv
  | n -> function Bounded ((lwb,upb),c) ->
          Bounded ((offset_bound n lwb, offset_bound n upb),.<.~c + n>.)
;;

(* Extract the dynamic value from bounded *)
let dyn_of_bounded = function
  | Bounded ((Low n,Low m),_) when n+1 = m -> .<n>.
  | Bounded (_,c) -> c
;;

let string_of_bounded : 'c bounded -> string = function
    Bounded (bounds,_) -> "B" ^ string_of_bounds bounds;;

(* Compare bounded: we compare bounds rather than code values *)
let same_bounded : 'c bounded -> 'c bounded -> bool =
  fun (Bounded (ab,_)) (Bounded (bb,_)) -> ab = bb;;

(* Check to see if one bounded immediately precedes the other.
   That is, the upper bound of one matches with the lower bound
   of the other.
*)
let successive_bounded : 'c bounded -> 'c bounded -> bool =
  fun (Bounded ((_,ab),_)) (Bounded ((bb,_),_)) -> ab = bb;;


let min_bounded = Bounded ((Low 0, Low 1),.<0>.);;
let array_bounds = (Low 0, High 0);;	(* Full bounds of an array *)

(* Check to see if a given bounded value is within the given bounds. 
   Return:
     - None if it is definitely out of bounds
     - Some c where c is the code (dynamic value) if it is definitely
       within the bounds
     - raise an error if cannot decide
*)
let within_bounds : 'c bounded -> bounds -> ('c,int) code option =
  fun (Bounded ((vl,vh),c)) (lwb, upb) ->
    if leq_bound vh lwb = True || leq_bound upb vl = True then None
    else if leq_bound lwb vl = True && leq_bound vh upb = True then Some c
    else failwith ("Can't decide if " ^ string_of_bounds (vl,vh) ^
		   " is within " ^ string_of_bounds (lwb,upb))
;;

(* A sample wrapper for our problems, which checks the upper
   bounds and generates the corresponding bounded value 
*)
let make_fun gen =
  .<fun n a b -> assert (n >= min_high); 
    .~(gen (Bounded ((High 0, High 1),.<n>.)) .<a>. .<b>.)>.
;;


(* Generate the for-loop *)
(* The loop variable i satisfies lwb <= i < upb *)
let forloop' lwb upb body =
  .<for i = .~lwb to (.~upb-1) do .~(body .<i>.) done>.;;

(* The same, with bounded index bounds *)
let forloop : 'c bounded -> 'c bounded ->
       ('c bounded -> ('c,unit) code) -> ('c,unit) code =
  fun (Bounded ((lwb,_),_) as lv) 
      (Bounded ((upb,_),_) as uv) body ->
  forloop' (dyn_of_bounded lv) (dyn_of_bounded uv) 
	  (fun j -> body (Bounded ((lwb,upb),j)))
;;
