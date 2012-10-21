(*
 Stencil computations: Takayuki Muranushi's challenge

 Signal processing or finite-element algorithms can often be stated
 as element-wise computations on shifted arrays. The following
 is the running example in the challenge.

  w = a * S[1] a
  b = a - w + S[-1] w

 Here |a| is a global input array, |b| is a global output array,
 and |w| is a working array.  The operation |S[n]|
 shifts the argument vector by |n| (to the left, if |n| is
 positive). All arithmetic operations on vectors (addition,
 multiplication, subtraction) are assumed elementwise. Global arrays
 are shared or distributed throughout a supercomputer; reading or writing
 them requires accessing off-the-chip memory or the inter-node
 communication.

 We should produce a program with a good or the optimal B/F (bytes-to-FLOPS
 ratio).  The key optimization is to cache the window of the input
 array (the stencil) and slide it. The problem is closely related to
 scalar promotion.

 To improve the B/F ratio, it helps to `do two time steps at once',
 (or, manually iterate the computation twice):

  w1 = a * S[1] a
  wm = a - w1 + S[-1] w1
  w2 = wm * S[1] wm
  b  = wm - w2 + S[-1] w2

 This file uses MetaOCaml with the delimcc library. Many parts of the
 present code are general and can be separated into a library.  The
 present code gradually develops a simple stencil DSL from scratch.
 See bounded.ml and bounded.mli for the library of partially known 
 dynamic values for safe array indexing.

The point of this code is to demonstrate how much one can develop
in a day (in Toukyou).
*)

open Bounded;;
(* 
  #load "bounded.cmo";;
*)

(* Configuration data *)

let filler = 0.;;	(* The filler for the values out of the array: *)
			(* halo values *)


(*
  Representation of a (read-only) array
  In general, an array is a function int -> float, or
  in our case, int -> float code.
  The index expression has the form of a statically known 
  offset from a dynamic value. We don't know the precise value
  of the dynamic value, but we do know its bounds.
  We will call such arrays `garr' (generalized arrays)
*)

type 'c index = I of int * 'c bounded;;
   (* 'c is an environment classifier *)

(* Offset the index *)
let plus : 'c index -> int -> 'c index = 
  fun i j -> match i with
    | I (x,c) -> I (x+j, c)
;;

(* Convert the index to the code *)
let idyn : 'c index -> ('c,int) code = function I (x,bv) ->
  dyn_of_bounded (offset_bounded x bv)
;;

let string_of_index : 'c index -> string = function 
  | I (x,bv) ->
      "I " ^ string_of_int x ^ "+" ^ string_of_bounded bv
;;

type 'c garr = 'c index -> ('c,float) code
;;

(* Convert a regular array to garr *)
let gref0 a : 'c garr = function I (x,bv) ->
  let bv = offset_bounded x bv in
  match within_bounds bv array_bounds with
  | None   -> .<filler>.
  | Some c -> .< (.~a).(.~c) >.;;

let gset0 a v = function I (x,bv) as ix ->
  let bv = offset_bounded x bv in
  match within_bounds bv array_bounds with
  | None   -> failwith ("output-of-bound ref to the output array: " ^
			string_of_index ix)
  | Some c -> .< (.~a).(.~c) <- .~v>.;;

(* Element-wise operations on garr: syntactic sugar *)
let ( +@ ) : 'c garr -> 'c garr -> 'c garr = fun a b j -> 
  .<.~(a j) +. .~(b j)>.;;
let ( -@ ) : 'c garr -> 'c garr -> 'c garr = fun a b j -> 
  .<.~(a j) -. .~(b j)>.;;
let ( *@ ) : 'c garr -> 'c garr -> 'c garr = fun a b j -> 
  .<.~(a j) *. .~(b j)>.;;

let ( <-@ ) : ('c, float array) code -> 'c garr -> 
              ('c index -> ('c,unit) code) = fun a b j -> 
  gset0 a (b j) j;;

(* The shift operation by a statically known amount *)
let ashift : int -> 'c garr -> 'c garr = fun s a -> function
  | I (x,b) -> a (I (x+s,b))
;;

  
(* Generate for-loop with halos: unroll the first halo and the last
   halo iterations.
   We re-use the previously written generator forloop.
   Keep in mind the right-to-left evaluation order of OCaml!
*)
let forlooph halo lwb upb body =
  let rec loop base j n = 		(* full loop unrolling *)
    if j >= n then .<()>. else
    let ix = I (0,offset_bounded j base) in
    if j = n-1 then body ix
    else let c1 = body ix in
	 let c2 = loop base (j+1) n in 
	 .<begin .~c1; .~c2 end>.
  in
  let c1 = loop lwb 0 halo in
  let lwb' = offset_bounded halo lwb and
      upb' = offset_bounded (-halo) upb in
  let c2 = forloop lwb' upb' (fun ib -> (body (I (0,ib)))) in
  let c3 = loop upb' 0 halo in
  .<begin .~c1; .~c2; .~c3 end>.
;;

(* Step I.
   The most naive code for the array computation, using halo of 1
*)

let t1 max_bounded a b =
  let a = gref0 a in
  forlooph 1 min_bounded max_bounded (
    let w = a *@ ashift 1 a in
    b <-@ a -@ w +@ ashift (-1) w)
;;

let t1c = make_fun t1;;


(* Generated code. It can be easily mapped to C or other low-level
   language.
   The code does handle boundary cases well; all array indices
   are in-bounds.
   The code is very unoptimal, with many, repeated input array
   references (e.g., a_2.(i_4) is repeated thrice).
*)

(*
val t1c : ('a, int -> float array -> float array -> unit) code =
  .<fun n_1 ->
     fun a_2 ->
      fun b_3 ->
       assert (n_1 >= 4);
       b_3.(0) <- ((a_2.(0) -. (a_2.(0) *. a_2.(0 + 1))) +. (0. *. a_2.(0)));
       for i_4 = 1 to ((n_1 + (-1)) - 1) do
        b_3.(i_4) <-
         ((a_2.(i_4) -. (a_2.(i_4) *. a_2.(i_4 + 1))) +.
           (a_2.(i_4 + (-1)) *. a_2.(i_4)))
       done;
       b_3.(n_1 + (-1)) <-
        ((a_2.(n_1 + (-1)) -. (a_2.(n_1 + (-1)) *. 0.)) +.
          (a_2.((n_1 + (-1)) + (-1)) *. a_2.(n_1 + (-1))))>.
*)


(* sample arrays *)
let a1 = Array.init 10 (fun i -> 1. +. 0.1 *. float_of_int i);;
let t1r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  (.! t1c) (Array.length a1) a1 b1; b1
;;
(* Check the result *)
let true = t1r =
  [|-0.100000000000000089; 0.880000000000000115; 0.96; 1.04000000000000026;
    1.12; 1.19999999999999929; 1.2799999999999998; 1.36000000000000032;
    1.44000000000000061; 5.32|]
;;

(* In contrast, if we don't use the halo, we get a problem *)
let t10 max_bounded a b =
  let a = gref0 a in
  forlooph 0 min_bounded max_bounded (
    let w = a *@ ashift 1 a in
    b <-@ a -@ w +@ ashift (-1) w)
;;

(*
let t10c = make_fun t10;;
Exception: Failure "Can't decide if [L-1,H-1) is within [L,H)".
*)

(* If the generator knows for sure that the index is out of bounds of the
   input array, we return the filler value. If the index is certain
   to be in-bounds, we generate the array reference.
   Without the halo, the index could be out-of-bounds (for the first and
   the last iteration of the loop), and could be in-bounds.
   The generator complains. It could've inserted the array-bound check.
   But we rather not add code to (inner) loops. 
   Unrolling the first and the last iteration makes the array access
   pattern certain, without additional checks.
*)

(* Step II.
   Now we develop a more sophisticated version, caching the
   accesses to main memory arrays.
   The approach is almost the same as optimizing Gibonacci,
   as described in our JFP paper.
*)

open Delimcc
;;

(* Memoize garr *)
let memo p (a : 'c garr) : 'c garr  =
 let table = ref [] in
 let deref j = 
   try List.assoc j !table 
     with Not_found ->
       shift p (fun k -> 
	 .<let t = .~(a j) in 
	 .~(table := (j,.<t>.) :: !table; k .<t>.) >.)
 in deref
;;

(* Add a prompt around the body of the loop.
   The prompt indicates the loop `scope' 
*)
let with_body_prompt p body = fun j -> 
  push_prompt p (fun () -> body j)
;;

  
(* The user code essentially remains the same. 
   We add the memoization, within the loop scope, as indicated by p.
   The prompt is created before forlooph runs its body; but the
   prompt is pushed only when forlooph runs its body.
*)
let t2 max_bounded a b =
  let p = new_prompt () in
  forlooph 1 min_bounded max_bounded (with_body_prompt p (
        let a = memo p (gref0 a) in
	let w = memo p (a *@ ashift 1 a) in
	b <-@ a -@ w +@ ashift (-1) w))
;;

let t2c = make_fun t2;;

(* The generated code: a_11.(i_18) array reference appears only once:
   we load the element into a local scalar t_20 and use the scalar.
   We still have to access the neighboring elements from the main
   memory. Therefore, B/F is very bad.
   We'll fixed that next.
*)
(*
val t2c : ('a, int -> float array -> float array -> unit) code =
  .<fun n_10 ->
     fun a_11 ->
      fun b_12 ->
       assert (n_10 >= 4);
       let t_14 = a_11.(0) in
       let t_15 = 0. in
       let t_13 = (t_15 *. t_14) in
       let t_17 = a_11.(0 + 1) in
       let t_16 = (t_14 *. t_17) in b_12.(0) <- ((t_14 -. t_16) +. t_13);
       for i_18 = 1 to ((n_10 + (-1)) - 1) do
        let t_20 = a_11.(i_18) in
        let t_21 = a_11.(i_18 + (-1)) in
        let t_19 = (t_21 *. t_20) in
        let t_23 = a_11.(i_18 + 1) in
        let t_22 = (t_20 *. t_23) in b_12.(i_18) <- ((t_20 -. t_22) +. t_19)
       done;
       let t_25 = a_11.(n_10 + (-1)) in
       let t_26 = a_11.((n_10 + (-1)) + (-1)) in
       let t_24 = (t_26 *. t_25) in
       let t_28 = 0. in
       let t_27 = (t_25 *. t_28) in
       b_12.(n_10 + (-1)) <- ((t_25 -. t_27) +. t_24)>.
*)

let t2r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  (.! t2c) (Array.length a1) a1 b1; b1
;;
let true = t1r = t2r;;
(* The sample result is the same *)


(*
  Step III.
  Stencil as a cache
*)

type 'c stencil = 
  'c bounded *				(* The base of the index i *)
  (int * ('c,float ref) code) list	(* Cache for a[i+j] for a range of j.
					   The list is ordered in increasing j.
					*)
;;

(* Insert a new association in order *)
let rec insert_in_order : 
  'k -> 'v -> ('k * 'v) list -> ('k * 'v) list = fun k v -> function 
  | []                      -> [(k,v)]
  | (k',v')::t when k <= k' -> (k,v)::(k',v')::t
  | kv :: t                 -> kv :: insert_in_order k v t
;;

(* Generate the code to slide the stencil, using garr to fill
   in the vacant slot. iv is already shifted.
*)

let stencil_slide : 'c stencil -> 'c garr -> ('c,unit) code = 
  fun (iv,table) a ->
  let rec loop = function
  | [(j,v)] -> .<.~v := .~(a (I (j,iv)))>.
  | ((_,vp)::(j,v)::rest) -> 
      .<begin .~vp := !(.~v); .~(loop ((j,v)::rest)) end>.
  in loop table
;;


(* stencil_memo transforms a garr into a memoized
  garr that uses a stencil as a cache. The signature of stencil_memo
  is essentially the same as that of memo.
  The function takes two prompts: prompt p is for the scope of a loop,
  and prompt q denotes a scope outside of the loop. References are allocated
  outside the loop but are assigned to inside.
  We have to watch for the evaluation order! 
  In OCaml bytecode, it is righ-to-left!
*)

let stencil_memo q p (a_base : 'c garr) : 'c garr  =
  let stencil = ref None in		(* uninitialized *)
  let make_ref ix =			(* outside the loop *)
    shift q (fun k ->			(* scope extrusion may happen! *)
      let ar = a_base ix in
      .<let s = ref .~ar in .~(k .<s>.)>.) in
					(* Cache access  *)
  let aref (inited, (iv',table)) ((I (j,iv)) as ix) = 
    try List.assoc j table 
    with Not_found -> 
     if inited then failwith ("A new index not seen at initialization: " ^ 
                              string_of_index ix)
     else  let ar = make_ref ix in
           stencil := Some (false, (iv',insert_in_order j ar table));
           ar
  in					(* dispatcher *)
  let check ((I (j,iv)) as ix) = 
    match !stencil with
    | None ->			(* initialization phase *)
	let ar = make_ref ix in
        stencil := Some (false,(iv,[(j,ar)])); ar
    | Some ((inited,(iv',table)) as st) when same_bounded iv' iv -> aref st ix
    | Some (inited,(iv',table)) when successive_bounded iv' iv ->
       shift p (fun k ->
         let sv = stencil_slide (iv,table) a_base in
         .<begin .~sv; .~(k ()) end>.);
       let st = (true,(iv,table)) in (* Freeze the stencil *)
       stencil := Some st; aref st ix
    | Some (inited,(iv',_))  ->
	failwith ("break in the loop? base  " ^ string_of_bounded iv' ^
                  " index " ^ string_of_index ix)
  in
  fun ix -> .<! .~(check ix)>.
;;

(* A loop that uses stencil to cache the accesses to the main array *)

let t3 max_bounded a b =
  let p = new_prompt () and q = new_prompt () in
  push_prompt q (fun () ->
  forlooph 1 min_bounded max_bounded (with_body_prompt p (
        let a = stencil_memo q p (gref0 a) in
	let w = (* stencil_memo q p *) (a *@ ashift 1 a) in
	b <-@ a -@ w +@ ashift (-1) w)))
;;

let t3c = make_fun t3;;

(* The code uses only one input array access per iteration,
   but the number of operations is a bit high.
   We should also use the stencil to cache the w computations.
   We fix this next.
*)
(*
val t3c : ('a, int -> float array -> float array -> unit) code =
  .<fun n_30 ->
     fun a_31 ->
      fun b_32 ->
       assert (n_30 >= 4);
       let s_33 = (ref a_31.(0)) in
       let s_34 = (ref 0.) in
       let s_35 = (ref a_31.(0 + 1)) in
       b_32.(0) <-
        (((! s_33) -. ((! s_33) *. (! s_35))) +. ((! s_34) *. (! s_33)));
       for i_36 = 1 to ((n_30 + (-1)) - 1) do
        begin
         (s_34 := (! s_33));
         (s_33 := (! s_35));
         (s_35 := a_31.(i_36 + 1))
        end;
        b_32.(i_36) <-
         (((! s_33) -. ((! s_33) *. (! s_35))) +. ((! s_34) *. (! s_33)))
       done;
       begin (s_34 := (! s_33)); (s_33 := (! s_35)); (s_35 := 0.) end;
       b_32.(n_30 + (-1)) <-
        (((! s_33) -. ((! s_33) *. (! s_35))) +. ((! s_34) *. (! s_33)))>.
*)

let t3r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  (.! t3c) (Array.length a1) a1 b1; b1
;;
let true = t1r = t3r;;
(* The sample result is the same *)

(* We can cache both the real array a and the virtual array 'w' *)
let t31 max_bounded a b =
  let p = new_prompt () and q = new_prompt () in
  push_prompt q (fun () ->
  forlooph 1 min_bounded max_bounded (with_body_prompt p (
        let a = stencil_memo q p (gref0 a) in
	let w = stencil_memo q p (a *@ ashift 1 a) in
	b <-@ a -@ w +@ ashift (-1) w)))
;;

let t31c = make_fun t31;;

(*
val t31c : ('a, int -> float array -> float array -> unit) code =
  .<fun n_37 ->
     fun a_38 ->
      fun b_39 ->
       assert (n_37 >= 4);
       let s_40 = (ref a_38.(0)) in
       let s_41 = (ref 0.) in
       let s_42 = (ref ((! s_41) *. (! s_40))) in
       let s_43 = (ref a_38.(0 + 1)) in
       let s_44 = (ref ((! s_40) *. (! s_43))) in
       b_39.(0) <- (((! s_40) -. (! s_44)) +. (! s_42));
       for i_45 = 1 to ((n_37 + (-1)) - 1) do
        begin
         (s_41 := (! s_40));
         (s_40 := (! s_43));
         (s_43 := a_38.(i_45 + 1))
        end;
        begin (s_42 := (! s_44)); (s_44 := ((! s_40) *. (! s_43))) end;
        b_39.(i_45) <- (((! s_40) -. (! s_44)) +. (! s_42))
       done;
       begin (s_41 := (! s_40)); (s_40 := (! s_43)); (s_43 := 0.) end;
       begin (s_42 := (! s_44)); (s_44 := ((! s_40) *. (! s_43))) end;
       b_39.(n_37 + (-1)) <- (((! s_40) -. (! s_44)) +. (! s_42))>.
*)

let t31r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  (.! t31c) (Array.length a1) a1 b1; b1
;;
let true = t1r = t31r;;
(* The sample result is the same *)



(* To improve B/F, we merge the consecutive simulation loops.
  We use a different stencil computation.

  w1 = a * S[1] a
  wm = a - w1 + S[-1] w1
  w2 = wm * S[1] wm
  b  = wm - w2 + S[-1] w2
*)

(* Do two time steps with the old code *)
let t20r = 
  let n = Array.length a1 in
  let b1 = Array.make n (-100.0) in
  (.! t1c) n a1 b1; 
  let c1 = Array.make n (-100.0) in
  (.! t1c) n b1 c1;
c1
;;
let true = t20r =
  [|-0.0119999999999999968; -0.052800000000000083; 0.806399999999999895;
    0.873599999999999932; 0.940800000000001191; 1.00799999999999979;
    1.0751999999999986; 1.14239999999999919; -4.26240000000000219;
    12.9808000000000039|]
;;

(* The baseline, the most naive fused code *)
(* The halo value must be two; otherwise, the generator complains *)
let t21 max_bounded a b =
  forlooph 2 min_bounded max_bounded (
    let a = gref0 a in
    let w1 = a *@ ashift 1 a in 
    let wm = a -@ w1 +@ ashift (-1) w1 in
    let w2 = wm *@ ashift 1 wm in
    b <-@ wm -@ w2 +@ ashift (-1) w2)
;;

let t21c = make_fun t21;;

(*
val t21c : ('a, int -> float array -> float array -> unit) code =
  .<fun n_46 ->
     fun a_47 ->
      fun b_48 ->
       assert (n_46 >= 4);
       begin
        b_48.(0) <-
         ((((a_47.(0) -. (a_47.(0) *. a_47.(0 + 1))) +. (0. *. a_47.(0))) -.
            (((a_47.(0) -. (a_47.(0) *. a_47.(0 + 1))) +. (0. *. a_47.(0)))
              *.
              ((a_47.(0 + 1) -. (a_47.(0 + 1) *. a_47.(0 + 2))) +.
                (a_47.(0) *. a_47.(0 + 1))))) +.
           (((0. -. (0. *. a_47.(0))) +. (0. *. 0.)) *.
             ((a_47.(0) -. (a_47.(0) *. a_47.(0 + 1))) +. (0. *. a_47.(0)))));
        b_48.(0 + 1) <-
         ((((a_47.(0 + 1) -. (a_47.(0 + 1) *. a_47.((0 + 1) + 1))) +.
             (a_47.((0 + 1) + (-1)) *. a_47.(0 + 1))) -.
            (((a_47.(0 + 1) -. (a_47.(0 + 1) *. a_47.((0 + 1) + 1))) +.
               (a_47.((0 + 1) + (-1)) *. a_47.(0 + 1))) *.
              ((a_47.((0 + 1) + 1) -.
                 (a_47.((0 + 1) + 1) *. a_47.((0 + 1) + 2))) +.
                (a_47.(0 + 1) *. a_47.((0 + 1) + 1))))) +.
           (((a_47.((0 + 1) + (-1)) -.
               (a_47.((0 + 1) + (-1)) *. a_47.(0 + 1))) +.
              (0. *. a_47.((0 + 1) + (-1)))) *.
             ((a_47.(0 + 1) -. (a_47.(0 + 1) *. a_47.((0 + 1) + 1))) +.
               (a_47.((0 + 1) + (-1)) *. a_47.(0 + 1)))))
       end;
       for i_49 = 2 to ((n_46 + (-2)) - 1) do
        b_48.(i_49) <-
         ((((a_47.(i_49) -. (a_47.(i_49) *. a_47.(i_49 + 1))) +.
             (a_47.(i_49 + (-1)) *. a_47.(i_49))) -.
            (((a_47.(i_49) -. (a_47.(i_49) *. a_47.(i_49 + 1))) +.
               (a_47.(i_49 + (-1)) *. a_47.(i_49))) *.
              ((a_47.(i_49 + 1) -. (a_47.(i_49 + 1) *. a_47.(i_49 + 2))) +.
                (a_47.(i_49) *. a_47.(i_49 + 1))))) +.
           (((a_47.(i_49 + (-1)) -. (a_47.(i_49 + (-1)) *. a_47.(i_49))) +.
              (a_47.(i_49 + (-2)) *. a_47.(i_49 + (-1)))) *.
             ((a_47.(i_49) -. (a_47.(i_49) *. a_47.(i_49 + 1))) +.
               (a_47.(i_49 + (-1)) *. a_47.(i_49)))))
       done;
       b_48.(n_46 + (-2)) <-
        ((((a_47.(n_46 + (-2)) -.
             (a_47.(n_46 + (-2)) *. a_47.((n_46 + (-2)) + 1))) +.
            (a_47.((n_46 + (-2)) + (-1)) *. a_47.(n_46 + (-2)))) -.
           (((a_47.(n_46 + (-2)) -.
               (a_47.(n_46 + (-2)) *. a_47.((n_46 + (-2)) + 1))) +.
              (a_47.((n_46 + (-2)) + (-1)) *. a_47.(n_46 + (-2)))) *.
             ((a_47.((n_46 + (-2)) + 1) -. (a_47.((n_46 + (-2)) + 1) *. 0.))
               +. (a_47.(n_46 + (-2)) *. a_47.((n_46 + (-2)) + 1))))) +.
          (((a_47.((n_46 + (-2)) + (-1)) -.
              (a_47.((n_46 + (-2)) + (-1)) *. a_47.(n_46 + (-2)))) +.
             (a_47.((n_46 + (-2)) + (-2)) *. a_47.((n_46 + (-2)) + (-1)))) *.
            ((a_47.(n_46 + (-2)) -.
               (a_47.(n_46 + (-2)) *. a_47.((n_46 + (-2)) + 1))) +.
              (a_47.((n_46 + (-2)) + (-1)) *. a_47.(n_46 + (-2))))));
       b_48.((n_46 + (-2)) + 1) <-
        ((((a_47.((n_46 + (-2)) + 1) -. (a_47.((n_46 + (-2)) + 1) *. 0.)) +.
            (a_47.(((n_46 + (-2)) + 1) + (-1)) *. a_47.((n_46 + (-2)) + 1)))
           -.
           (((a_47.((n_46 + (-2)) + 1) -. (a_47.((n_46 + (-2)) + 1) *. 0.))
              +.
              (a_47.(((n_46 + (-2)) + 1) + (-1)) *. a_47.((n_46 + (-2)) + 1)))
             *. ((0. -. (0. *. 0.)) +. (a_47.((n_46 + (-2)) + 1) *. 0.)))) +.
          (((a_47.(((n_46 + (-2)) + 1) + (-1)) -.
              (a_47.(((n_46 + (-2)) + 1) + (-1)) *. a_47.((n_46 + (-2)) + 1)))
             +.
             (a_47.(((n_46 + (-2)) + 1) + (-2)) *.
               a_47.(((n_46 + (-2)) + 1) + (-1)))) *.
            ((a_47.((n_46 + (-2)) + 1) -. (a_47.((n_46 + (-2)) + 1) *. 0.))
              +.
              (a_47.(((n_46 + (-2)) + 1) + (-1)) *. a_47.((n_46 + (-2)) + 1)))))>.
*)

(* Try the code with two fused time steps *)
let t21r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  (.! t21c) (Array.length a1) a1 b1; b1
;;
let true = t20r = t21r;;

(* Now we optimize the code with stencils *)
(* We merely add the stencil annotations *)
let t25 max_bounded a b =
  let p = new_prompt () and q = new_prompt () in
  push_prompt q (fun () ->
  forlooph 2 min_bounded max_bounded (with_body_prompt p (
    let a  = stencil_memo q p (gref0 a) in
    let w1 = stencil_memo q p (a *@ ashift 1 a) in 
    let wm = stencil_memo q p (a -@ w1 +@ ashift (-1) w1) in
    let w2 = stencil_memo q p (wm *@ ashift 1 wm) in
    b <-@ wm -@ w2 +@ ashift (-1) w2)))
;;

let t25c = make_fun t25;;



(* The code below has one read access and one write access per
   iteration, and does 6 FLOPS. So, B/F is 16/6, the optimal.
*)

(*
val t25c : ('a, int -> float array -> float array -> unit) code =
  .<fun n_50 ->
     fun a_51 ->
      fun b_52 ->
       assert (n_50 >= 4);
       let s_53 = (ref a_51.(0)) in
       let s_54 = (ref 0.) in
       let s_55 = (ref ((! s_54) *. (! s_53))) in
       let s_56 = (ref a_51.(0 + 1)) in
       let s_57 = (ref ((! s_53) *. (! s_56))) in
       let s_58 = (ref (((! s_53) -. (! s_57)) +. (! s_55))) in
       let s_59 = (ref 0.) in
       let s_60 = (ref ((! s_59) *. (! s_54))) in
       let s_61 = (ref (((! s_54) -. (! s_55)) +. (! s_60))) in
       let s_62 = (ref ((! s_61) *. (! s_58))) in
       let s_63 = (ref a_51.(0 + 2)) in
       let s_64 = (ref ((! s_56) *. (! s_63))) in
       let s_65 = (ref (((! s_56) -. (! s_64)) +. (! s_57))) in
       let s_66 = (ref ((! s_58) *. (! s_65))) in
       begin
        b_52.(0) <- (((! s_58) -. (! s_66)) +. (! s_62));
        begin
         (s_59 := (! s_54));
         (s_54 := (! s_53));
         (s_53 := (! s_56));
         (s_56 := (! s_63));
         (s_63 := a_51.((0 + 1) + 2))
        end;
        begin
         (s_60 := (! s_55));
         (s_55 := (! s_57));
         (s_57 := (! s_64));
         (s_64 := ((! s_56) *. (! s_63)))
        end;
        begin
         (s_61 := (! s_58));
         (s_58 := (! s_65));
         (s_65 := (((! s_56) -. (! s_64)) +. (! s_57)))
        end;
        begin (s_62 := (! s_66)); (s_66 := ((! s_58) *. (! s_65))) end;
        b_52.(0 + 1) <- (((! s_58) -. (! s_66)) +. (! s_62))
       end;
       for i_67 = 2 to ((n_50 + (-2)) - 1) do
        begin
         (s_59 := (! s_54));
         (s_54 := (! s_53));
         (s_53 := (! s_56));
         (s_56 := (! s_63));
         (s_63 := a_51.(i_67 + 2))
        end;
        begin
         (s_60 := (! s_55));
         (s_55 := (! s_57));
         (s_57 := (! s_64));
         (s_64 := ((! s_56) *. (! s_63)))
        end;
        begin
         (s_61 := (! s_58));
         (s_58 := (! s_65));
         (s_65 := (((! s_56) -. (! s_64)) +. (! s_57)))
        end;
        begin (s_62 := (! s_66)); (s_66 := ((! s_58) *. (! s_65))) end;
        b_52.(i_67) <- (((! s_58) -. (! s_66)) +. (! s_62))
       done;
       begin
        begin
         (s_59 := (! s_54));
         (s_54 := (! s_53));
         (s_53 := (! s_56));
         (s_56 := (! s_63));
         (s_63 := 0.)
        end;
        begin
         (s_60 := (! s_55));
         (s_55 := (! s_57));
         (s_57 := (! s_64));
         (s_64 := ((! s_56) *. (! s_63)))
        end;
        begin
         (s_61 := (! s_58));
         (s_58 := (! s_65));
         (s_65 := (((! s_56) -. (! s_64)) +. (! s_57)))
        end;
        begin (s_62 := (! s_66)); (s_66 := ((! s_58) *. (! s_65))) end;
        b_52.(n_50 + (-2)) <- (((! s_58) -. (! s_66)) +. (! s_62))
       end;
       begin
        (s_59 := (! s_54));
        (s_54 := (! s_53));
        (s_53 := (! s_56));
        (s_56 := (! s_63));
        (s_63 := 0.)
       end;
       begin
        (s_60 := (! s_55));
        (s_55 := (! s_57));
        (s_57 := (! s_64));
        (s_64 := ((! s_56) *. (! s_63)))
       end;
       begin
        (s_61 := (! s_58));
        (s_58 := (! s_65));
        (s_65 := (((! s_56) -. (! s_64)) +. (! s_57)))
       end;
       begin (s_62 := (! s_66)); (s_66 := ((! s_58) *. (! s_65))) end;
       b_52.((n_50 + (-2)) + 1) <- (((! s_58) -. (! s_66)) +. (! s_62))>.
*)

(* Try the code with two fused time steps *)
let t25r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  (.! t25c) (Array.length a1) a1 b1; b1
;;
let true = t20r = t25r;;


(*
Hand-written code
  w1 = a * S[1] a
  wm = a - w1 + S[-1] w1
  w2 = wm * S[1] wm
  b  = wm - w2 + S[-1] w2

*)

let stencil_hand n a b =
  let a_m2   = ref 0.0 in
  let a_m1   = ref 0.0 in
  let a_0    = ref a.(0) in
  let a_1    = ref a.(1) in
  let a_2    = ref a.(2) in
  let w1_m2  = ref (!a_m2 *. !a_m1) in
  let w1_m1  = ref (!a_m1 *. !a_0) in
  let w1_0   = ref (!a_0  *. !a_1) in
  let w1_1   = ref (!a_1  *. !a_2) in
  let wm_m1  = ref (!a_m1 -. !w1_m1 +. !w1_m2) in
  let wm_0   = ref (!a_0 -. !w1_0 +. !w1_m1) in
  let wm_1   = ref (!a_1 -. !w1_1 +. !w1_0) in
  let w2_m1  = ref (!wm_m1 *. !wm_0) in
  let w2_0   = ref (!wm_0  *. !wm_1) in
  b.(0) <- !wm_0 -. !w2_0 +. !w2_m1;

  a_1 := !a_2; a_2 := a.(3);
  w1_0 := !w1_1;  w1_1 := (!a_1  *. !a_2);
  wm_0 := !wm_1;  wm_1 := (!a_1 -. !w1_1 +. !w1_0);
  w2_m1 := !w2_0; w2_0 := (!wm_0  *. !wm_1);
  b.(1) <- !wm_0 -. !w2_0 +. !w2_m1;

  for i = 2 to (n-3) do
   a_1 := !a_2; a_2 := a.(i+2);
   w1_0 := !w1_1;  w1_1 := (!a_1  *. !a_2);
   wm_0 := !wm_1;  wm_1 := (!a_1 -. !w1_1 +. !w1_0);
   w2_m1 := !w2_0; w2_0 := (!wm_0  *. !wm_1);
   b.(i) <- !wm_0 -. !w2_0 +. !w2_m1;
  done;

  a_1 := !a_2; a_2 := 0.0;
  w1_0 := !w1_1;  w1_1 := (!a_1  *. !a_2);
  wm_0 := !wm_1;  wm_1 := (!a_1 -. !w1_1 +. !w1_0);
  w2_m1 := !w2_0; w2_0 := (!wm_0  *. !wm_1);
  b.(n-2) <- !wm_0 -. !w2_0 +. !w2_m1;
  a_1 := !a_2; a_2 := 0.0;
  w1_0 := !w1_1;  w1_1 := (!a_1  *. !a_2);
  wm_0 := !wm_1;  wm_1 := (!a_1 -. !w1_1 +. !w1_0);
  w2_m1 := !w2_0; w2_0 := (!wm_0  *. !wm_1);
  b.(n-1) <- !wm_0 -. !w2_0 +. !w2_m1;
;;


let stencil_hand_r = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  stencil_hand (Array.length a1) a1 b1; b1
;;

let true = t20r = stencil_hand_r;;


print_endline "\nAll done\n";;

(*
  Old code, kept for reference


(* We introduce a different format for array computation, 
   somewhat reminiscent of conjunctive normal form.
*)

(* A term is a constant factor and a sequence
   a_shift factors. Each factor is represented by its shift value.
   Thus a term is a constant and a list of shift factors.
   We keep the factors sorted.
*)

type term = Term of float * int list;;
type formula = term list;;

(* The running example
   dref 0 -@ dref 0 *@ dref 1 +@ dref (-1) *@ dref 0))
*)

let f1 = [
  Term (1., [0]);
  Term (-1., [0;1]);
  Term (1.,[-1;0])]
;;

let term_shift n (Term (c,shifts)) =
  Term (c,List.map (fun x -> x + n) shifts)
;;
let formula_shift n : formula -> formula
  = List.map (term_shift n);;

(* multiply two formulas, the result is unnormalized *)
let formula_mul f1 f2 : formula =
  let ft_mul (Term (c1,s1)) = 
   List.map (function Term (c2,s2) -> Term (c1*.c2,List.sort compare (s1 @ s2)))
  in
  List.concat (List.map (fun t1 -> ft_mul t1 f2) f1)
;;

let rec formulas_mul = function
  | []  -> []
  | [x] -> x
  | h1::h2::t -> formulas_mul (formula_mul h1 h2 :: t)
;;

(* substitute a formula for a term. The result is unnormalized:
   it may contain repeated terms
*)
let formula_subst_unnorm formula (Term (c,shifts)) : formula =
  List.map (function Term (c',s) -> Term (c'*.c,s)) (
   formulas_mul (List.map (fun n -> formula_shift n formula) shifts))
;;

(* Merge similar terms *)
let rec formula_normalize f : formula =
  let check = function
  | []  -> []
  | [t] -> [t]
  | Term (c1,s1)::Term (c2,s2)::rest when s1 = s2 ->
    let c = c1 +. c2 in
    if c = 0. then formula_normalize rest 
    else formula_normalize (Term (c,s1)::rest)
  | t1::rest -> t1 :: formula_normalize rest
  in 
  let cmp (Term (c1,s1)) (Term (c2,s2)) = compare s1 s2 in
  check (List.sort cmp f)
;;

let formula_conv f1 f2 : formula = 
  formula_normalize 
   (List.concat (List.map (fun t -> formula_subst_unnorm f1 t) f2))
;;

let f2 = formula_conv f1 f1;;

*)
