/*
 * The Yices SMT Solver. Copyright 2015 SRI International.
 *
 * This program may only be used subject to the noncommercial end user
 * license agreement which is downloadable along with this program.
 */
 
#pragma once

#include <stdio.h>

#ifndef NDEBUG
#define DD_DEBUG
#define DD_STATS
#endif

#include <cudd.h>

#include "terms/terms.h"
#include "utils/pointer_vectors.h"

typedef DdNode BDD;

#define BDDS_RESERVE_MAX 2

typedef struct bdd_manager_s bdd_manager_t;

/** Reference for BDD vectors. */
typedef uint32_t bddvec_id_t;

/** Null vector ID */
extern const bddvec_id_t bddvec_null_id;

/** BDD vector managed by the BDD vector manager */
typedef struct bddvec_s {
  uint32_t size;
  bddvec_id_t id;
} bddvec_t;

/** Actual null vector */
extern const bddvec_t bddvec_null;

/** Swap the two BDDs */
static inline
void bddvec_swap(bddvec_t* x, bddvec_t* y) {
  bddvec_t tmp = *x; *x = *y; *y = tmp;
}

/**
 * The structure responsible for bare BDD computation. It relies on the CUDD
 * library to perform the BDD computation.
 *
 * All BDD manipulating functions use the BDD manager to allocate the result.
 * This means that the BDD manager's BDD allocations need to be invariant.
 */
typedef struct {
  DdManager* cudd;
  struct bdd_manager_s* bddm;
  int* tmp_inputs;
  char* tmp_model;
  size_t tmp_alloc_size;
} CUDD;

/** Construct and allocate cudd (passed bddm will not be used in construction) */
CUDD* bdds_new(bdd_manager_t* bddm);

/** Destruct and delete cudd */
void bdds_delete(CUDD* cudd);

/**
 * Given the term and BDDs of all the children compute the BDDs into
 * the output.
 *
 * The vector children_bdds contains a pointer for each child. The pointer
 * is of the type BDD** and is of the same size as the term children.
 */
bddvec_t bdds_compute_bdds(CUDD* cudd, term_table_t* terms, term_t t, const pvector_t* children_bdds);

/** Initialize: set all to NULL. */
void bdds_init(BDD** a, uint32_t n);

/** Dereference all non-NULL bdds in a and set them to NULL */
void bdds_clear(CUDD* cudd, BDD** a, uint32_t n);

/** Dereference all non-NULL bdds  */
void bdds_detach(CUDD* cudd, BDD** a, uint32_t n);

/** Attach extra reference all bdds in a. */
void bdds_attach(BDD** a, uint32_t n);

/** Compare the two BDD vectors. */
bool bdds_eq(BDD** a, BDD** b, uint32_t n);

/** Print the BDDs to out. */
void bdds_print(CUDD* cudd, BDD** a, uint32_t n, FILE* out);

/** Check if the disjoint BDDs are a point of given size (only one solution). */
bool bdds_is_point(CUDD* cudd, BDD** a, uint32_t n, uint32_t bitsize);

/**
 * Check if the constant satisfies the constraint C(x). The variables in x
 * should be the same size as x_value.
 */
bool bdds_is_model(CUDD* cudd, BDD** x, BDD* C_x, const bvconstant_t* x_value);

/**
 * Get a constant that satisfies the constraint C(x). Constant will be changed
 * only in places that matter (i.e., you should init it with bit-vealues that
 * you want as default.
 */
void bdds_get_model(CUDD* cudd, BDD** x, BDD* C_x, bvconstant_t* out);

/** Make a new variable. */
void bdds_mk_variable(CUDD* cudd, BDD** out, uint32_t n);

/** Make a repeat b...b BDD. */
void bdds_mk_repeat(CUDD* cudd, BDD** out, BDD* b, uint32_t n);

/** Make a constant 0...0 BDD. */
void bdds_mk_zero(CUDD* cudd, BDD** out, uint32_t n);

  /** Make a constant 0...01 BDD. */
void bdds_mk_one(CUDD* cudd, BDD** out, uint32_t n);

/** Make a constant BDD. */
void bdds_mk_constant(CUDD* cudd, BDD** out, uint32_t n, const bvconstant_t* c);

/** Check if a BDD is constant */
bool bdds_is_constant(CUDD* cudd, BDD** a, uint32_t n);

/** Check if a BDD is constant 0...0 */
bool bdds_is_constant_zero(CUDD* cudd, BDD** a, uint32_t n);

/** Check if a BDD is constant 0...01 */
bool bdds_is_constant_one(CUDD* cudd, BDD** a, uint32_t n);

/** Check if a BDD is constant 1...1 */
bool bdds_is_constant_neg_one(CUDD* cudd, BDD** a, uint32_t n);

/** Check if a BDD is a power of 2. Returns power, or -1 if not */
int32_t bdds_is_constant_pow2(CUDD* cudd, BDD** a, uint32_t n);

/** Check whether the two BDDs are disjoint */
bool bdds_are_disjoint(CUDD* cudd, BDD* a, BDD* b);

/** Negate the BDDs a. */
void bdds_mk_not(CUDD* cudd, BDD** out, BDD** a, uint32_t n);

/** Boolean and of the BDDs in a and b. */
void bdds_mk_and(CUDD* cudd, BDD** out, BDD** a, BDD** b, uint32_t n);

/** Boolean and of the bdds in a, i.e. out[0] = a[0] && ... && a[n-1] */
void bdds_mk_conjunction(CUDD* cudd, BDD** out, BDD** a, uint32_t n);

/** Boolean or of the BDDs in a and b. */
void bdds_mk_or(CUDD* cudd, BDD** out, BDD** a, BDD** b, uint32_t n);

/** Two's complement (negation) of the BDDs in a. */
void bdds_mk_2s_complement(CUDD* cudd, BDD** out, BDD** a, uint32_t n);

/** Addition of the BDDs in a and b. */
void bdds_mk_plus(CUDD* cudd, BDD** out, BDD** a, BDD** b, uint32_t n);

/** Multiplication of the BDDs in a and b. */
void bdds_mk_mult(CUDD* cudd, BDD** out, BDD** a, BDD** b, uint32_t n);

/** Division of the BDDs in a and b. */
void bdds_mk_div(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Remainder of the BDDs in a and b. */
void bdds_mk_rem(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Signed division of BDDs in a and b. */
void bdds_mk_sdiv(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Signed remainder of BDDs in a and b (rounding to 0). */
void bdds_mk_srem(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Signed remainder of BDDs in a and b. (rounding to -infinity). */
void bdds_mk_smod(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Left shift of BDDs in a and b (padding 0). */
void bdds_mk_shl(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Left shift of BDDs in a */
void bdds_mk_shl_const(CUDD* cudd, BDD** out_bdds, BDD** a, uint32_t shift, uint32_t n);

/** Logical shift right of BDDs in a and b (padding with 0). */
void bdds_mk_lshr(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Arithmetic shift right (padding with sign bit). */
void bdds_mk_ashr(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Equality circuit of BDDs in a and b, out is of n 1. */
void bdds_mk_eq(CUDD* cudd, BDD** out, BDD** a, BDD** b, uint32_t n);

/** Unsigned comparison circuit of BDDs in a and b, out is of n 1. */
void bdds_mk_ge(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Signed comparison circuit of BDDs in a and b, out is of n 1. */
void bdds_mk_sge(CUDD* cudd, BDD** out_bdds, BDD** a, BDD** b, uint32_t n);

/** Add BDDs in to_add to a conjunction of disjoint BDDs. */
void bdds_disjoint_set_add(CUDD* cudd, BDD** a, uint32_t n, pvector_t* disjoint_set);
