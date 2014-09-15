(set-logic QF_AUFLIA)
(set-info :source |
  These benchmarks come from SMT-queries of the infinite-state
  model-checking tool SAFARI (see SAFARI: SMT-based Abstraction
  For Arrays with Interpolants, by F. Alberti, R. Bruttomesso,
  S. Ghilardi, S. Ranise, and N. Sharygina, in CAV 2012).
  Each benchmark is a full trace resulting from the verification
  of a program manipulating arrays (e.g., find and element in an
  array, sort an array, etc.)
  Generated by Roberto Bruttomesso.
|)
(set-info :smt-lib-version 2.0)
(set-info :category "industrial")
; Declaring all variables
(declare-fun x () Int)
; success
(declare-fun j () Int)
; success
(declare-fun z0 () Int)
; success
(declare-fun z1 () Int)
; success
(declare-fun z2 () Int)
; success
(declare-fun z3 () Int)
; success
(declare-fun z4 () Int)
; success
(declare-fun z5 () Int)
; success
(declare-fun z6 () Int)
; success
; All INDEX variables are non-negative
(assert (<= 0 z0) )
; success
(assert (<= 0 z1) )
; success
(assert (<= 0 z2) )
; success
(assert (<= 0 z3) )
; success
(assert (<= 0 z4) )
; success
(assert (<= 0 z5) )
; success
(assert (<= 0 z6) )
; success
; All variables are different
(assert (not (= z0 z1)) )
; success
(assert (not (= z0 z2)) )
; success
(assert (not (= z0 z3)) )
; success
(assert (not (= z0 z4)) )
; success
(assert (not (= z0 z5)) )
; success
(assert (not (= z0 z6)) )
; success
(assert (not (= z1 z2)) )
; success
(assert (not (= z1 z3)) )
; success
(assert (not (= z1 z4)) )
; success
(assert (not (= z1 z5)) )
; success
(assert (not (= z1 z6)) )
; success
(assert (not (= z2 z3)) )
; success
(assert (not (= z2 z4)) )
; success
(assert (not (= z2 z5)) )
; success
(assert (not (= z2 z6)) )
; success
(assert (not (= z3 z4)) )
; success
(assert (not (= z3 z5)) )
; success
(assert (not (= z3 z6)) )
; success
(assert (not (= z4 z5)) )
; success
(assert (not (= z4 z6)) )
; success
(assert (not (= z5 z6)) )
; success
; Declaring non INDEX variables
; Declaring scalar variables
(declare-fun I () Int)
; success
(declare-fun J () Int)
; success
(declare-fun f () Int)
; success
(declare-fun l () Int)
; success
(declare-fun a_length () Int)
; success
; Declaring array variables
(declare-fun a () (Array Int Int))
; success
; Asserting 'globality' of some arrays
; Asserting system axiom
(assert (and (<= 1 l) (<= l 2)) )
; success
; Asserting system axiom
(assert (not (<= a_length 0)) )
; success
(set-info :status sat)
(check-sat)
; sat
(push 1)
; success
(define-fun Initial () Bool (and (= I 0) (= J 0) (= f 0) (= l 1)))
; success
(define-fun Node0_0 ((z0 Int)) Bool (and (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (= l 2)) )
; success
(push 1)
; success
; 
; 
; 
; Check initial intersection
(push 1)
; success
(assert (and (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (= l 2)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (= l 2)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (= l 2)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; End fixpoint check
; 
; 
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (= l 2)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; End fixpoint check
; 
; 
(push 1)
; success
(assert (and (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (= l 2)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; Checking if the new preimage computed from node 0 with t1(z0) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (= z0 I)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; Checking if the new preimage computed from node 0 with t1(z1) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; 
; Defining new node derived from node 1 applying transition 1
(define-fun Node1_0 ((z0 Int) (z1 Int)) Bool (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; End fixpoint check
; 
; 
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; End fixpoint check
; 
; 
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_0 () Int)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (<= ta_0 I) (not (<= ta_0 z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (<= ta_0 I) (not (<= ta_0 z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_1 () Int)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (<= ta_1 I) (not (<= ta_1 z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (<= ta_1 I) (not (<= ta_1 z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_4 () Int)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_5 () Int)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_6 () Int)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (not (<= ta_6 z0)) (= z1 ta_6)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_7 () Int)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (not (<= ta_7 z0)) (= z1 ta_7)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_8 () Int)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_9 () Int)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_10 () Int)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0)) (= ta_10 1)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0)) (= ta_10 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_11 () Int)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0)) (= ta_11 1)) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0)) (= ta_11 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_12 () Int)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; BEGIN TERM ABSTRACTION
(declare-fun ta_13 () Int)
; success
; END TERM ABSTRACTION
(push 1)
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; 
; Defining new node derived from node 1 applying transition 1
(define-fun Node1_1 ((z0 Int) (z1 Int)) Bool (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
; 
; 
; Checking valid implication
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (<= a_length z0)) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I)) )
; success
(assert (not (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) ) )
; success
(assert (not (and (<= 0 z1) (not (= z0 z1)) (= z0 I) (not (= (select a z1) 0)) (not (<= I z1))) ) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; 
; 
; 
; 
; 
; Check initial intersection
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (and (= I 0) (= J 0) (= f 0) (= l 1)) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (Node0_0 z0 )) )
; success
(assert (not (Node0_0 z1 )) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; End fixpoint check
; 
; 
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (Node0_0 z0 )) )
; success
(assert (not (Node0_0 z1 )) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; End fixpoint check
; 
; 
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; Checking if the new preimage computed from node 1 with t0(z0) is safe
; Checking if the new preimage computed from node 1 with t0(z1) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (<= a_length I)) (not (= z0 z1)) (= z1 I) (not (<= (+ I 1) z0)) (= z1 (+ I 1))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; Checking if the new preimage computed from node 1 with t0(z2) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (<= a_length I)) (not (= z1 z2)) (not (= z0 z1)) (not (= z0 z2)) (= z2 I) (not (<= (+ I 1) z0)) (= z1 (+ I 1))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
(push 1)
; success
(assert (and (<= 0 z0) (not (= (select a z0) 0)) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; Checking if the new preimage computed from node 1 with t1(z0) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z0 I) (= z1 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; Checking if the new preimage computed from node 1 with t1(z1) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; Checking if the new preimage computed from node 1 with t1(z2) is safe
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (<= a_length I) (not (= z1 z2)) (not (= z0 z1)) (not (= z0 z2)) (= z1 I) (= z2 I) (not (<= I z0))) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; 
; Defining new node derived from node 2 applying transition 0
(define-fun Node2_0 ((z0 Int) (z1 Int) (z2 Int)) Bool (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (<= a_length I)) (not (= z1 z2)) (not (= z0 z1)) (not (= z0 z2)) (= z2 I) (not (<= (+ I 1) z0)) (= z1 (+ I 1))) )
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (<= a_length I)) (not (= z1 z2)) (not (= z0 z1)) (not (= z0 z2)) (= z2 I) (not (<= (+ I 1) z0)) (= z1 (+ I 1))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (not (<= a_length I)) (not (= z1 z2)) (not (= z0 z1)) (not (= z0 z2)) (= z2 I) (not (<= (+ I 1) z0)) (= z1 (+ I 1))) )
; success
(assert (not (Node1_1 z0 z1 )) )
; success
(assert (not (Node1_1 z0 z2 )) )
; success
(assert (not (Node1_1 z1 z0 )) )
; success
(assert (not (Node1_1 z1 z2 )) )
; success
(assert (not (Node1_1 z2 z0 )) )
; success
(assert (not (Node1_1 z2 z1 )) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; End fixpoint check
; 
; Defining new node derived from node 3 applying transition 1
(define-fun Node3_0 ((z0 Int) (z1 Int)) Bool (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(push 1)
; success
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(set-info :status sat)
(check-sat)
; sat
(pop 1)
; success
; 
; 
; Fixpoint test starts (3)
(push 1)
; success
; Asserting the label of the node
(assert (and (= l 1) (<= 0 z0) (not (= (select a z0) 0)) (<= a_length I) (not (= z0 z1)) (= z1 I) (not (<= I z0))) )
; success
(assert (not (Node1_1 z0 z1 )) )
; success
(assert (not (Node1_1 z1 z0 )) )
; success
(set-info :status unsat)
(check-sat)
; unsat
(pop 1)
; success
; End fixpoint check
(pop 1)
; success
(pop 1)
; success
(exit)
