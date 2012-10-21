(* Hand-written reference code for the Stencil challenge and the unit tests *)

(* The code is based on the problem posed by Takayuki Muranushi.
   The code below properly handles the edge cases.
*)

(* Specification

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

*)

let filler = 0.;;	(* The filler for the values out of the array: *)
			(* halo values *)

let spec_1step : int -> float array -> float array -> unit = fun n a b ->
  for i=0 to n-1 do
    let a j = if j >= 0 && j < n then a.(j) else filler in
    let w j = a j *. a (j+1) in
    b.(i) <- a i -. w i +. w (i-1)
done
;;

(* sample arrays *)
let a1 = Array.init 10 (fun i -> 1. +. 0.1 *. float_of_int i);;
let t1 = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  spec_1step (Array.length a1) a1 b1; b1
;;

(* Check the result *)
let true = t1 =
  [|-0.100000000000000089; 0.880000000000000115; 0.96; 1.04000000000000026;
    1.12; 1.19999999999999929; 1.2799999999999998; 1.36000000000000032;
    1.44000000000000061; 5.32|]
;;

(* Repeat spec_1step twice: apply spec_1step to the input array a
   producing b. Then apply spec_1step to b producing c.
*)

let t20 = 
  let n = Array.length a1 in
  let b1 = Array.make n (-100.0) in
  spec_1step n a1 b1; 
  let c1 = Array.make n (-100.0) in
  spec_1step n b1 c1;
c1
;;
(* Unit test *)
let true = t20 =
  [|-0.0119999999999999968; -0.052800000000000083; 0.806399999999999895;
    0.873599999999999932; 0.940800000000001191; 1.00799999999999979;
    1.0751999999999986; 1.14239999999999919; -4.26240000000000219;
    12.9808000000000039|]
;;

(*
Hand-written optimized code for the two-step computation

  w1 = a * S[1] a
  wm = a - w1 + S[-1] w1
  w2 = wm * S[1] wm
  b  = wm - w2 + S[-1] w2

The challenge is to produce the code automatically, by applying
the stencil optimization to the above specification.
*)

let optimized_2step n a b =
  let a_m2   = ref filler in
  let a_m1   = ref filler in
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

  a_1 := !a_2; a_2 := filler;
  w1_0 := !w1_1;  w1_1 := (!a_1  *. !a_2);
  wm_0 := !wm_1;  wm_1 := (!a_1 -. !w1_1 +. !w1_0);
  w2_m1 := !w2_0; w2_0 := (!wm_0  *. !wm_1);
  b.(n-2) <- !wm_0 -. !w2_0 +. !w2_m1;
  a_1 := !a_2; a_2 := filler;
  w1_0 := !w1_1;  w1_1 := (!a_1  *. !a_2);
  wm_0 := !wm_1;  wm_1 := (!a_1 -. !w1_1 +. !w1_0);
  w2_m1 := !w2_0; w2_0 := (!wm_0  *. !wm_1);
  b.(n-1) <- !wm_0 -. !w2_0 +. !w2_m1;
;;


let t25 = 
  let b1 = Array.make (Array.length a1) (-100.0) in
  optimized_2step (Array.length a1) a1 b1; b1
;;

let true = t20 = t25;;


print_endline "\nAll done\n";;
