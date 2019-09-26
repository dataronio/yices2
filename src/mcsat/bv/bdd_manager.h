/*
 * This file is part of the Yices SMT Solver.
 * Copyright (C) 2017 SRI International.
 *
 * Yices is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Yices is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Yices.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <stdio.h>

#include "terms/terms.h"
#include "terms/bv_constants.h"
#include "mcsat/mcsat_types.h"

#include "bdd_computation.h"

/**
 * Structure responsible for all BDD interactions.
 *
 * Intended use:
 * - create the manager;
 * - register all bit-vector variables with the manager;
 * - create bdd's for terms that are unit over added variables
 * - compute with bdd's (e.g., and)
 * - pick bit-vector value from a BDD
 * - allocate/deallocate new temporary BDDs that you need
 */
typedef struct bdd_manager_s bdd_manager_t;

/** Create a new BDD manager (allocate and construct) */
bdd_manager_t* bdd_manager_new(const plugin_context_t* ctx);

/** Delete the given BDD manager (destruct and deallocate) */
void bdd_manager_delete(bdd_manager_t* bddm);

/** Constant function true */
bddvec_t bdd_manager_true(const bdd_manager_t* bddm);

/** Constant function false */
bddvec_t bdd_manager_false(const bdd_manager_t* bddm);

/** Add term to the BDD manager (only added terms can be set to a value) */
void bdd_manager_add_term(bdd_manager_t* bddm, term_t t);

/** Check if BDD manager is aware of the term (has been added) */
bool bdd_manager_has_term(const bdd_manager_t* bddm, term_t term);

/** Set a value for a variable (copies the value) */
void bdd_manager_set_bv_value(bdd_manager_t* bddm, term_t x, const bvconstant_t* value);

/** Set a value for a 1-bit variable */
void bdd_manager_set_bool_value(bdd_manager_t* bddm, term_t t, bool value);

/**
 * Get a BDD representing a given Boolean literal. All variables except for
 * x will be replaced with their set values.
 *
 * @param t a term of the form t(x, y1, ..., yn), where all terms y_i have set values
 * @param t x a variable whose value will not be used
 *
 * @return bdds a BDD vector (conjunction) representing the constraint [refcount attached]
 */
bddvec_t bdd_manager_get(bdd_manager_t* bddm, term_t bv_literal, term_t x);

/** Get the actuall BDDs of the vector (use with care, and changes to m migth reallocate) */
BDD** bdd_manager_get_bdds(const bdd_manager_t* bddm, bddvec_t v);

/** Pick a value for a given variable x that satisfies the given v = BDD(x). */
void bdd_manager_pick_value(bdd_manager_t* bddm, term_t x, bddvec_t v, bvconstant_t* out);

/** Check if a given value satisfies the given BDD(x) */
bool bdd_manager_is_model(bdd_manager_t* bddm, term_t x, bddvec_t v, const bvconstant_t* x_value);

/** Detach a BDD vector */
void bdd_manager_detach(bdd_manager_t* bddm, bddvec_t v);

/** Detach a BDD vector and set all to NULL */
void bdd_manager_delete_vec(bdd_manager_t* bddm, bddvec_t v);

/** Make a copy of the given vector, returns attached version */
bddvec_t bdd_manager_new_copy(bdd_manager_t* bddm, bddvec_t v);

/** Attach a BDD vector */
void bdd_manager_attach(bdd_manager_t* bddm, bddvec_t v);

/** Print BDDs to output */
void bdd_manager_print_bddvec(const bdd_manager_t* bddm, bddvec_t v, FILE* out);

/** Check whether the given BDD is an empty set */
bool bdd_manager_is_empty(const bdd_manager_t* bddm, bddvec_t v);

/** Check whether the given BDD represents one solution. */
bool bdd_manager_is_point(const bdd_manager_t* bddm, bddvec_t bdd, uint32_t bitsize);

/** Intersect the two BDDs (result attached) */
bddvec_t bdd_manager_intersect(bdd_manager_t* bddm, bddvec_t v1, bddvec_t v2);
