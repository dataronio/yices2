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

#include <mcsat/bv/bdd_manager.h>
#include "bdd_computation.h"
#include "bv_utils.h"
#include "utils/memalloc.h"
#include "utils/int_hash_map.h"
#include "utils/int_hash_sets.h"
#include "utils/ptr_vectors.h"
#include "utils/int_vectors.h"
#include "terms/terms.h"

#include "mcsat/plugin.h"
#include "mcsat/tracing.h"

#define BDDVEC_ALLOCATOR_VECTORS_PER_CHUNK 100

/** Manage fixed size allocations */
typedef struct bddvec_allocator_s {
  /** Size of the bitvector */
  uint32_t vec_size;
  /** Memory for the BDDs */
  pvector_t chunks;
  /** List of free id's */
  ivector_t free_list;
  /** Top id */
  bddvec_id_t next_to_allocate;
} bddvec_allocator_t;

void bddvec_allocator_construct(bddvec_allocator_t* m, uint32_t size) {
  m->vec_size = size;
  m->next_to_allocate = 0;
  init_pvector(&m->chunks, 0);
  init_ivector(&m->free_list, 0);
}

void bddvec_allocator_destruct(bddvec_allocator_t* m) {
  delete_ivector(&m->free_list);
  for (uint32_t i = 0; i < m->chunks.size; ++ i) {
    BDD** bdds_i = (BDD**) m->chunks.data[i];
    safe_free(bdds_i);
  }
  delete_pvector(&m->chunks);
}

BDD** bddvec_allocator_get(bddvec_allocator_t* m, bddvec_id_t id) {
  uint32_t chunk_i = id / BDDVEC_ALLOCATOR_VECTORS_PER_CHUNK;
  uint32_t i = id % BDDVEC_ALLOCATOR_VECTORS_PER_CHUNK;
  assert(chunk_i < m->chunks.size);
  BDD** bdds = ((BDD**) m->chunks.data[chunk_i]) + (i*m->vec_size);
  return bdds;
}

void bddvec_allocator_delete_vec(bddvec_allocator_t* m, bddvec_id_t id) {
  BDD** bdds = bddvec_allocator_get(m, id);
  for (BDD** i = bdds; i != bdds + m->vec_size; ++ i) {
    *i = NULL;
  }
  ivector_push(&m->free_list, id);
}

bddvec_id_t bddvec_allocator_new_vec(bddvec_allocator_t* m) {

  bddvec_id_t result = bddvec_null_id;

  // Check free list
  if (m->free_list.size > 0) {
    result = ivector_pop2(&m->free_list);
  } else {
    // New id is at the end
    result = m->next_to_allocate ++;
    // Check if we need a new chunk
    uint32_t chunk_i = result / BDDVEC_ALLOCATOR_VECTORS_PER_CHUNK;
    if (chunk_i == m->chunks.size) {
      uint32_t bdds_in_chunk = m->vec_size * BDDVEC_ALLOCATOR_VECTORS_PER_CHUNK;
      BDD** new_chunk = safe_malloc(sizeof(BDD*) * bdds_in_chunk);
      for (BDD** bdd = new_chunk; bdd != new_chunk + bdds_in_chunk; ++ bdd) {
        *bdd = NULL;
      }
      pvector_push(&m->chunks, new_chunk);
    }
  }
  return result;
}

typedef struct bddvec_manager_s {
  uint32_t size;
  uint32_t capacity;
  bddvec_allocator_t* alloc;
} bddvec_manager_t;

#define bddvec_MANAGER_DEFAULT_SIZE 64;

void bddvec_manager_ensure_allocated(bddvec_manager_t* m) {
  uint32_t i;
  uint32_t n = m->capacity;
  while (m->size > n) {
    n += 1;
    n += n >> 2;
  }
  if (n > m->capacity) {
    m->alloc = (bddvec_allocator_t*) safe_realloc(m->alloc, n*sizeof(bddvec_allocator_t));
    for (i = m->capacity; i < n; ++ i) {
      bddvec_allocator_construct(m->alloc + i, i);
    }
    m->capacity = n;
  }
}

void bddvec_manager_ensure_size(bddvec_manager_t* m, uint32_t size) {
  if (size >= m->size) {
    m->size = size + 1;
    bddvec_manager_ensure_allocated(m);
  }
}

void bddvec_manager_construct(bddvec_manager_t* m) {
  m->alloc = NULL;
  m->size = bddvec_MANAGER_DEFAULT_SIZE;
  m->capacity = 0;
  bddvec_manager_ensure_allocated(m);
}

void bddvec_manager_destruct(bddvec_manager_t* m) {
  uint32_t i;
  for (i = 0; i < m->capacity; ++ i) {
    bddvec_allocator_destruct(m->alloc + i);
  }
  safe_free(m->alloc);
}

bddvec_t bddvec_manager_new_vec(bddvec_manager_t* m, uint32_t size) {
  bddvec_manager_ensure_size(m, size);
  bddvec_t result = { size, bddvec_allocator_new_vec(m->alloc + size) };
  return result;
}

BDD** bddvec_manager_get_bdds(const bddvec_manager_t* m, bddvec_t v) {
  assert(v.size < m->size);
  return bddvec_allocator_get(m->alloc + v.size, v.id);
}

void bddvec_manager_delete_vec(bddvec_manager_t* m, bddvec_t v) {
  assert(v.size < m->size);
  bddvec_allocator_delete_vec(m->alloc + v.size, v.id);
}

#define BDDVEC_NULL_ID UINT32_MAX

/** Null vector id */
const bddvec_id_t bddvec_null_id = BDDVEC_NULL_ID;

/** Actual null vector */
const bddvec_t bddvec_null = { 0, BDDVEC_NULL_ID };

/** Information stored about the terms. */
typedef struct {
  /** The term (size of the term is in value.size) */
  term_t t;
  /** Id of the BDDs (only valid if timestamp != 0) */
  bddvec_t v;
  /** If of the constant BDDs */
  bddvec_t v_const;
  /** Time-stamp of the last BDD change modulo unassigned variable (0 for never) */
  uint32_t bdd_timestamp;
  /** Unit variable of term t, when BDD was computed (NULL if none) */
  term_t unassigned_variable;
  /** The bit-vector value of the term */
  bvconstant_t value;
  /** Time-stamp of the last value change (0 for never) */
  uint32_t value_timestamp;
} term_info_t;

static
void term_info_construct(term_info_t* term_info, term_t t, bddvec_t v, bddvec_t v_const, uint32_t bitsize) {
  assert(bitsize > 0);
  term_info->t = t;
  term_info->v = v;
  term_info->v_const = v_const;
  term_info->bdd_timestamp = 0;
  term_info->unassigned_variable = NULL_TERM;
  init_bvconstant(&term_info->value);
  bvconstant_set_bitsize(&term_info->value, bitsize);
  term_info->value_timestamp = 0;
}

/** All data needed for managing the BDDs. */
struct bdd_manager_s {
  /** Term table */
  const plugin_context_t* ctx;
  /** Cudd manager to be used for all */
  CUDD* cudd;

  /** Map from terms to their index in the term_info array */
  int_hmap_t term_to_info_index;

  /** Info on the terms */
  term_info_t* term_info;
  /** Size of the info */
  uint32_t term_info_size;
  /** Capacity of the info */
  uint32_t term_info_capacity;
  /** Timestamp of the last value change */
  uint32_t timestamp;

  /** Vector of all the terms terms added */
  ivector_t term_list;

  /** Memory for storing the BDD pointers */
  bddvec_manager_t bdds;

  /** The variable that we consider unassigned when computing */
  term_t unassigned_var;

  /** Marks for visited nodes */
  int_hmap_t visited;

  /** List of terms to recompute the value */
  ivector_t value_recompute;

  /** List of terms to recompute the BDD */
  ivector_t bdd_recompute;

  /** BDD constant false */
  BDD* bdd_false;

  /** BDD constant true */
  BDD* bdd_true;

  /** BDD constant false in vector form */
  bddvec_t bdd_false_v;

  /** BDD constant true in vector form */
  bddvec_t bdd_true_v;

  /** Temporaray BDD storage */
  pvector_t bdd_temp;
};

/** Delete all BDD info (constant, and the BDDs allocated) */
static
void term_info_destruct(term_info_t* term_info) {
  delete_bvconstant(&term_info->value);
}

bdd_manager_t* bdd_manager_new(const plugin_context_t* ctx) {
  // Allocate
  bdd_manager_t* bddm = (bdd_manager_t*) safe_malloc(sizeof(bdd_manager_t));

  bddm->ctx= ctx;
  bddm->cudd = bdds_new(bddm);
  bddm->term_info = NULL;
  bddm->term_info_size = 0;
  bddm->term_info_capacity = 0;
  bddm->timestamp = 1;

  bddvec_manager_construct(&bddm->bdds);

  init_int_hmap(&bddm->term_to_info_index, 0);
  init_ivector(&bddm->term_list, 0);

  bddm->unassigned_var = NULL_TERM;

  init_int_hmap(&bddm->visited, 0);
  init_ivector(&bddm->value_recompute, 0);
  init_ivector(&bddm->bdd_recompute, 0);

  bddm->bdd_false = NULL;
  bddm->bdd_true = NULL;
  bdds_mk_zero(bddm->cudd, &bddm->bdd_false, 1);
  bdds_mk_one(bddm->cudd, &bddm->bdd_true, 1);

  bddm->bdd_false_v = bddvec_manager_new_vec(&bddm->bdds, 1);
  BDD** bdd_false_v_bdds = bddvec_manager_get_bdds(&bddm->bdds, bddm->bdd_false_v);
  bdd_false_v_bdds[0] = bddm->bdd_false;

  bddm->bdd_true_v = bddvec_manager_new_vec(&bddm->bdds, 1);
  BDD** bdd_true_v_bdds = bddvec_manager_get_bdds(&bddm->bdds, bddm->bdd_true_v);
  bdd_true_v_bdds[0] = bddm->bdd_true;

  init_pvector(&bddm->bdd_temp, 0);

  return bddm;
}

bddvec_t bdd_manager_true(const bdd_manager_t* bddm) {
  return bddm->bdd_true_v;
}

bddvec_t bdd_manager_false(const bdd_manager_t* bddm) {
  return bddm->bdd_false_v;
}

static inline
term_info_t* bdd_manager_get_info(const bdd_manager_t* bddm, term_t t) {
  // Destruct the term info
  int_hmap_pair_t* info_find = int_hmap_find((int_hmap_t*) &bddm->term_to_info_index, t);
  assert(info_find != NULL);
  uint32_t t_info_index = info_find->val;
  term_info_t* t_info = bddm->term_info + t_info_index;
  return t_info;
}

/**
 * If the term contains the unassigned variable: return term BDDs
 * Otherwise: return constant BDDs
 */
static inline
BDD** bdd_manager_get_bdds_same_size(bdd_manager_t* bddm, term_t t) {
  term_info_t* t_info = bdd_manager_get_info(bddm, t);
  int_hmap_pair_t* find = int_hmap_find(&bddm->visited, t);
  uint32_t t_bitsize = bv_term_bitsize(bddm->ctx->terms, t);
  assert(find != NULL);
  if (find->val) {
    // Contains the unassigned variable, so we get the BDDs
    if (t_info->v.size != t_bitsize) {
      assert(t_bitsize == 1);
      // Have to flatten the BDD representation and replace it
      bddvec_t old_v = t_info->v;
      bddvec_t new_v = bddvec_manager_new_vec(&bddm->bdds, 1);
      BDD** old_bdds = bddvec_manager_get_bdds(&bddm->bdds, old_v);
      BDD** new_bdds = bddvec_manager_get_bdds(&bddm->bdds, new_v);
      bdds_mk_conjunction(bddm->cudd, new_bdds, old_bdds, old_v.size);
      bddvec_manager_delete_vec(&bddm->bdds, old_v);
      t_info->v = new_v;
    }
    return bddvec_manager_get_bdds(&bddm->bdds, t_info->v);
  } else {
    // Convert the value to BDD
    const bvconstant_t* value = &t_info->value;
    BDD** out = bddvec_manager_get_bdds(&bddm->bdds, t_info->v_const);
    bdds_clear(bddm->cudd, out, value->bitsize);
    bdds_mk_constant(bddm->cudd, out, value->bitsize, value);
    return out;
  }
}

void bdd_manager_delete(bdd_manager_t* bddm) {

  CUDD* cudd = bddm->cudd;

  // Decrease reference counts and destruct all term info
  for (uint32_t i = 0; i < bddm->term_list.size; ++ i) {
    // Term and info
    term_t t = bddm->term_list.data[i];
    term_info_t* t_info = bdd_manager_get_info(bddm, t);
    bdd_manager_delete_vec(bddm, t_info->v);
    bddvec_manager_delete_vec(&bddm->bdds, t_info->v);
    bdd_manager_delete_vec(bddm, t_info->v_const);
    bddvec_manager_delete_vec(&bddm->bdds, t_info->v_const);
    term_info_destruct(t_info);
  }

  bddvec_manager_destruct(&bddm->bdds);

  delete_int_hmap(&bddm->term_to_info_index);
  delete_ivector(&bddm->term_list);

  delete_int_hmap(&bddm->visited);
  delete_ivector(&bddm->value_recompute);
  delete_ivector(&bddm->bdd_recompute);

  bdds_clear(cudd, &bddm->bdd_false, 1);
  bdds_clear(cudd, &bddm->bdd_true, 1);

  bdds_delete(cudd);

  delete_pvector(&bddm->bdd_temp);

  safe_free(bddm->term_info);
  safe_free(bddm);
}

#define BV_BDD_MEM_DEFAULT_SIZE 10
#define BV_BDD_MEM_MAX_SIZE (UINT32_MAX/8)

static
void bdd_manager_ensure_term_capacity(bdd_manager_t* bddm, uint32_t var_index) {
  uint32_t n;
  while (bddm->term_info_capacity <= var_index) {
    n = bddm->term_info_capacity;
    if (n == 0) {
      n = BV_BDD_MEM_DEFAULT_SIZE;
    } else {
      n ++;
      n += n >> 1;
      if (n >= BV_BDD_MEM_MAX_SIZE) {
        out_of_memory();
      }
    }
    bddm->term_info = (term_info_t*) safe_realloc(bddm->term_info, n*sizeof(term_info_t));
    bddm->term_info_capacity = n;
  }
}

bool bdd_manager_has_term(const bdd_manager_t* bddm, term_t term) {
  int_hmap_pair_t* find = int_hmap_find((int_hmap_t*) &bddm->term_to_info_index, term);
  return (find != NULL);
}

uint32_t bdd_manager_allocate_term_info(bdd_manager_t* bddm, term_t t, uint32_t bitsize) {

  assert(int_hmap_find(&bddm->term_to_info_index, t) == NULL);
  assert(bv_term_bitsize(bddm->ctx->terms, t) == bitsize);

  // Add the reference for t
  uint32_t t_info_index = bddm->term_info_size;
  int_hmap_add(&bddm->term_to_info_index, t, t_info_index);
  // Make sure there is enough memory
  bdd_manager_ensure_term_capacity(bddm, t_info_index);
  // Increase the size
  term_info_t* x_info = bddm->term_info + t_info_index;
  bddm->term_info_size ++;
  // Add the BDD vector for t
  bddvec_t v = bddvec_manager_new_vec(&bddm->bdds, bitsize);
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  bdds_init(v_bdds, v.size);
  // Add the BDD vector for t's constant value
  bddvec_t v_const = bddvec_manager_new_vec(&bddm->bdds, bitsize);
  BDD** v_const_bdds = bddvec_manager_get_bdds(&bddm->bdds, v_const);
  bdds_init(v_const_bdds, v_const.size);
  // Construct info
  term_info_construct(x_info, t, v, v_const, bitsize);
  // Remember the allocated terms
  ivector_push(&bddm->term_list, t);

  return t_info_index;
}

/**
 * Recursively allocate the data for all bit-vector terms in t.
 */
static
void bdd_manager_ensure_term_data(bdd_manager_t* bddm, term_t t, uint32_t bitsize) {

  term_table_t* terms = bddm->ctx->terms;

  int_hmap_pair_t* info_index_find = int_hmap_find(&bddm->term_to_info_index, t);
  if (info_index_find == NULL) {

    if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
      ctx_trace_printf(bddm->ctx, "bdd_manager: allocating term/BDD data for ");
      ctx_trace_term(bddm->ctx, t);
    }

    // Allocate this term
    uint32_t t_info_index = bdd_manager_allocate_term_info(bddm, t, bitsize);
    term_info_t* t_info = bddm->term_info + t_info_index;
    BDD** t_bdds = bddvec_manager_get_bdds(&bddm->bdds, t_info->v);

    // Ensure data for the sub-terms
    if (is_neg_term(t)) {
      term_t t_pos = unsigned_term(t);
      bdd_manager_ensure_term_data(bddm, t_pos, bitsize);
    } else {
      term_kind_t t_kind = term_kind(terms, t);
      switch (t_kind) {
      case OR_TERM: // Boolean
      case EQ_TERM: // Boolean equality
      case BV_EQ_ATOM:
      case BV_GE_ATOM:
      case BV_SGE_ATOM: {
        // Boolean atoms, 2 children are bitvectors
        assert(bitsize == 1);
        composite_term_t* atom_comp = composite_term_desc(terms, t);
        for (uint32_t i = 0; i < atom_comp->arity; ++ i) {
          uint32_t child_bitsize = bv_term_bitsize(terms, atom_comp->arg[i]);
          bdd_manager_ensure_term_data(bddm, atom_comp->arg[i], child_bitsize);
        }
        break;
      }
      case BV_ARRAY:
      case BV_DIV:
      case BV_REM:
      case BV_SDIV:
      case BV_SREM:
      case BV_SMOD:
      case BV_SHL:
      case BV_LSHR:
      case BV_ASHR: {
        // Ensure data for children
        composite_term_t* t_comp = composite_term_desc(terms, t);
        for (uint32_t i = 0; i < t_comp->arity; ++ i) {
          uint32_t t_i_bitsize = bv_term_bitsize(terms, t_comp->arg[i]);
          assert(t_kind == BV_ARRAY || t_i_bitsize == bitsize);
          assert(t_kind != BV_ARRAY || t_i_bitsize == 1);
          bdd_manager_ensure_term_data(bddm, t_comp->arg[i], t_i_bitsize);
        }
        break;
      }
      case BIT_TERM: {
        select_term_t* desc = bit_term_desc(terms, t);
        term_t child = desc->arg;
        uint32_t child_bitsize = bv_term_bitsize(terms, child);
        bdd_manager_ensure_term_data(bddm, child, child_bitsize);
        break;
      }
      case BV_POLY: {
        // Polynomial, allocate data for the subterms (variables)
        bvpoly_t* t_poly = bvpoly_term_desc(terms, t);
        for (uint32_t i = 0; i < t_poly->nterms; ++i) {
          if (t_poly->mono[i].var == const_idx) continue;
          term_t var = t_poly->mono[i].var;
          assert(bv_term_bitsize(terms, var) == bitsize);
          bdd_manager_ensure_term_data(bddm, var, bitsize);
        }
        break;
      }
      case BV64_POLY: {
        // Polynomial, allocate data for the subterms (variables)
        bvpoly64_t* t_poly = bvpoly64_term_desc(terms, t);
        for (uint32_t i = 0; i < t_poly->nterms; ++i) {
          if (t_poly->mono[i].var == const_idx) continue;
          term_t var = t_poly->mono[i].var;
          assert(bv_term_bitsize(terms, var) == bitsize);
          bdd_manager_ensure_term_data(bddm, var, bitsize);
        }
        break;
      }
      case POWER_PRODUCT: {
        pprod_t* t_pprod = pprod_term_desc(terms, t);
        for (uint32_t i = 0; i < t_pprod->len; ++ i) {
          term_t var = t_pprod->prod[i].var;
          assert(bv_term_bitsize(terms, var) == bitsize);
          bdd_manager_ensure_term_data(bddm, var, bitsize);
        }
        break;
      }
      case CONSTANT_TERM:
        if (t == true_term) {
          bvconst_set_bit(t_info->value.data, 0);
        } else if (t == false_term) {
          bvconst_clr_bit(t_info->value.data, 0);
        } else {
          assert(false); // Only Boolean constants
        }
        t_info->value_timestamp = 1;
        break;
      case BV_CONSTANT: {
        // Set the value
        bvconst_term_t* t_desc = bvconst_term_desc(terms, t);
        bvconstant_copy(&t_info->value, t_desc->bitsize, t_desc->data);
        t_info->value_timestamp = 1;
        break;
      }
      case BV64_CONSTANT: {
        // Set the value
        bvconst64_term_t* t_desc = bvconst64_term_desc(terms, t);
        bvconstant_copy64(&t_info->value, t_desc->bitsize, t_desc->value);
        t_info->value_timestamp = 1;
        break;
      }
      default:
        // We get here only with variables, or foreign terms.
        assert(bv_term_kind_get_type(t_kind) == BV_TERM_VARIABLE);
        // Add the variables nodes
        bdds_mk_variable(bddm->cudd, t_bdds, bitsize);
        // Mark the info for this variable
        t_info->unassigned_variable = t;
        t_info->bdd_timestamp = 1;
        break;
      }
    }
  } else {
    assert(info_index_find != NULL);
  }
}

void bdd_manager_add_term(bdd_manager_t* bddm, term_t t) {
  uint32_t bitsize = bv_term_bitsize(bddm->ctx->terms, t);
  bdd_manager_ensure_term_data(bddm, t, bitsize);
}

/**
 * Go through the term and:
 * - compute the value and value timestamp of each subterm that evaluates
 * - compute the BDD timestamp of each subterm
 * - if the BDD needs to be recomputed, put it in the list
 *
 * Returns true if the term contains the unassigned variable.
 *
 * If the return value is true the value timestamp should be ignored.
 * If the return value is false the BDD timestamp should be ignored.
 */
static
bool bdd_manager_recompute_timestamps(bdd_manager_t* bddm, term_t t, uint32_t bitsize, uint32_t* bdd_timestamp, uint32_t* value_timestamp) {

  term_table_t* terms = bddm->ctx->terms;

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: recomputing timestamps for ");
    ctx_trace_term(bddm->ctx, t);
  }

  // Info for t
  term_info_t* t_info = bdd_manager_get_info(bddm, t);

  // If visited already, skip
  int_hmap_pair_t* find = int_hmap_find(&bddm->visited, t);
  if (find != NULL) {
    *bdd_timestamp = t_info->bdd_timestamp;
    *value_timestamp = t_info->value_timestamp;
    return find->val;
  }

  // Does this term contain the unassigned variable
  bool contains_unassigned = false;

  // Variables: timestamps and values are always OK, we set them manually
  if (bv_term_is_variable(terms, t)) {
    *value_timestamp = t_info->value_timestamp;
    *bdd_timestamp = t_info->bdd_timestamp;
    assert(*bdd_timestamp == 1);
    contains_unassigned = (t == bddm->unassigned_var);
    int_hmap_add(&bddm->visited, t, contains_unassigned);
    return contains_unassigned;
  }

  // Initialize the timestamps
  *bdd_timestamp = 0;
  *value_timestamp = 0;

  // Should we recompute value
  bool recompute_value = false;

  // Should we recompute BDD
  bool recompute_bdd = false;

  // Negation
  if (is_neg_term(t)) {
    contains_unassigned = bdd_manager_recompute_timestamps(bddm, unsigned_term(t), bitsize, bdd_timestamp, value_timestamp);
  } else {
    term_kind_t t_kind = term_kind(terms, t);
    switch (t_kind) {
    case OR_TERM: // Boolean OR
    case EQ_TERM: // Boolean equality
    case BV_EQ_ATOM:
    case BV_GE_ATOM:
    case BV_SGE_ATOM:
    case BV_ARRAY:
    case BV_DIV:
    case BV_REM:
    case BV_SDIV:
    case BV_SREM:
    case BV_SMOD:
    case BV_SHL:
    case BV_LSHR:
    case BV_ASHR: {
      // Get the children BDDs
      composite_term_t* t_comp = composite_term_desc(terms, t);
      for (uint32_t i = 0; i < t_comp->arity; ++i) {
        uint32_t bdd_timestamp_i = 0;
        uint32_t value_timestamp_i = 0;
        uint32_t t_i = t_comp->arg[i];
        uint32_t bitsize_i = bv_term_bitsize(terms, t_i);
        bool t_i_contains_unassigned = bdd_manager_recompute_timestamps(bddm, t_i, bitsize_i, &bdd_timestamp_i, &value_timestamp_i);
        contains_unassigned = contains_unassigned || t_i_contains_unassigned;
        if (bdd_timestamp_i > *bdd_timestamp) { *bdd_timestamp = bdd_timestamp_i; }
        if (value_timestamp_i > *value_timestamp) { *value_timestamp = value_timestamp_i; }
      }
      break;
    }
    case BIT_TERM:
      contains_unassigned = bdd_manager_recompute_timestamps(bddm, bit_term_arg(terms, t), 1, bdd_timestamp, value_timestamp);
      break;
    case BV_POLY: {
      bvpoly_t* t_poly = bvpoly_term_desc(terms, t);
      for (uint32_t i = 0; i < t_poly->nterms; ++i) {
        if (t_poly->mono[i].var == const_idx) continue;
        uint32_t bdd_timestamp_i = 0;
        uint32_t value_timestamp_i = 0;
        term_t t_i = t_poly->mono[i].var;
        uint32_t bitsize_i = bv_term_bitsize(terms, t_i);
        bool t_i_contains_unassigned = bdd_manager_recompute_timestamps(bddm, t_i, bitsize_i, &bdd_timestamp_i, &value_timestamp_i);
        contains_unassigned = contains_unassigned || t_i_contains_unassigned;
        if (bdd_timestamp_i > *bdd_timestamp) { *bdd_timestamp = bdd_timestamp_i; }
        if (value_timestamp_i > *value_timestamp) { *value_timestamp = value_timestamp_i; }
      }
      break;
    }
    case BV64_POLY: {
      // Polynomial, allocate data for the subterms (variables)
      bvpoly64_t* t_poly = bvpoly64_term_desc(terms, t);
      for (uint32_t i = 0; i < t_poly->nterms; ++i) {
        if (t_poly->mono[i].var == const_idx) continue;
        uint32_t bdd_timestamp_i = 0;
        uint32_t value_timestamp_i = 0;
        term_t t_i = t_poly->mono[i].var;
        uint32_t bitsize_i = bv_term_bitsize(terms, t_i);
        bool t_i_contains_unassigned = bdd_manager_recompute_timestamps(bddm, t_i, bitsize_i, &bdd_timestamp_i, &value_timestamp_i);
        contains_unassigned = contains_unassigned || t_i_contains_unassigned;
        if (bdd_timestamp_i > *bdd_timestamp) { *bdd_timestamp = bdd_timestamp_i; }
        if (value_timestamp_i > *value_timestamp) { *value_timestamp = value_timestamp_i; }
      }
      break;
    }
    case POWER_PRODUCT: {
      pprod_t* t_pprod = pprod_term_desc(terms, t);
      for (uint32_t i = 0; i < t_pprod->len; ++ i) {
        uint32_t bdd_timestamp_i = 0;
        uint32_t value_timestamp_i = 0;
        term_t t_i = t_pprod->prod[i].var;
        uint32_t bitsize_i = bv_term_bitsize(terms, t_i);
        bool t_i_contains_unassigned = bdd_manager_recompute_timestamps(bddm, t_i, bitsize_i, &bdd_timestamp_i, &value_timestamp_i);
        contains_unassigned = contains_unassigned || t_i_contains_unassigned;
        if (bdd_timestamp_i > *bdd_timestamp) { *bdd_timestamp = bdd_timestamp_i; }
        if (value_timestamp_i > *value_timestamp) { *value_timestamp = value_timestamp_i; }
      }
      break;
    }
    case CONSTANT_TERM:
    case BV_CONSTANT:
    case BV64_CONSTANT:
      // Nothing to do, always constant
      *bdd_timestamp = 1;
      *value_timestamp = 1;
      contains_unassigned = false;
      break;
    default:
      // Shouldn't be here, variables handled above
      assert(false);
    }
  }

  // Should we recompute
  recompute_value = !contains_unassigned &&
      (t_info->value_timestamp != *value_timestamp);
  recompute_bdd = contains_unassigned &&
      (t_info->bdd_timestamp != *bdd_timestamp || t_info->unassigned_variable != bddm->unassigned_var);
  // Set the timestamps
  if (recompute_value) {
    t_info->value_timestamp = *value_timestamp;
    ivector_push(&bddm->value_recompute, t);
  }
  if (recompute_bdd) {
    t_info->bdd_timestamp = *bdd_timestamp;
    ivector_push(&bddm->bdd_recompute, t);
  }

  // Mark as visited
  int_hmap_add(&bddm->visited, t, contains_unassigned);

  return contains_unassigned;
}

static inline
bvconstant_t* bdd_manager_get_value(bdd_manager_t* bddm, term_t t) {
  term_info_t* t_info = bdd_manager_get_info(bddm, t);
  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "value of ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), t_info->value.data, t_info->value.bitsize);
    ctx_trace_printf(bddm->ctx, "\n");
  }
  return &t_info->value;
}

static inline
void bdd_manager_compute_value(bdd_manager_t* bddm, term_t t) {

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: computing value for ");
    ctx_trace_term(bddm->ctx, t);
  }

  uint32_t i;
  term_t t_i;
  bvconstant_t* value_i;
  term_table_t* terms = bddm->ctx->terms;

  pvector_t children_values;
  init_pvector(&children_values, 0);

  // Negation
  if (is_neg_term(t)) {
    t_i = unsigned_term(t);
    value_i = bdd_manager_get_value(bddm, t_i);
    pvector_push(&children_values, value_i);
  } else {
    term_kind_t t_kind = term_kind(terms, t);
    switch (t_kind) {
    case OR_TERM:
    case EQ_TERM:
    case BV_EQ_ATOM:
    case BV_GE_ATOM:
    case BV_SGE_ATOM:
    case BV_ARRAY:
    case BV_DIV:
    case BV_REM:
    case BV_SDIV:
    case BV_SREM:
    case BV_SMOD:
    case BV_SHL:
    case BV_LSHR:
    case BV_ASHR: {
      // Get the children BDDs
      composite_term_t* t_comp = composite_term_desc(terms, t);
      for (i = 0; i < t_comp->arity; ++i) {
        t_i = t_comp->arg[i];
        value_i = bdd_manager_get_value(bddm, t_i);
        pvector_push(&children_values, value_i);
      }
      break;
    }
    case BIT_TERM:
      t_i = bit_term_arg(terms, t);
      value_i = bdd_manager_get_value(bddm, t_i);
      pvector_push(&children_values, value_i);
      break;
    case BV_POLY: {
      bvpoly_t* t_poly = bvpoly_term_desc(terms, t);
      for (uint32_t i = 0; i < t_poly->nterms; ++i) {
        if (t_poly->mono[i].var == const_idx) continue;
        t_i = t_poly->mono[i].var;
        value_i = bdd_manager_get_value(bddm, t_i);
        pvector_push(&children_values, value_i);
      }
      break;
    }
    case BV64_POLY: {
      // Polynomial, allocate data for the subterms (variables)
      bvpoly64_t* t_poly = bvpoly64_term_desc(terms, t);
      for (uint32_t i = 0; i < t_poly->nterms; ++i) {
        if (t_poly->mono[i].var == const_idx) continue;
        t_i = t_poly->mono[i].var;
        value_i = bdd_manager_get_value(bddm, t_i);
        pvector_push(&children_values, value_i);
      }
      break;
    }
    case POWER_PRODUCT:{
      pprod_t* t_pprod = pprod_term_desc(terms, t);
      for (uint32_t i = 0; i < t_pprod->len; ++ i) {
        term_t t_i = t_pprod->prod[i].var;
        value_i = bdd_manager_get_value(bddm, t_i);
        pvector_push(&children_values, value_i);
      }
      break;
    }
    default:
      // Shouldn't be here
      assert(false);
    }
  }

  // We have the children values compute
  term_info_t* t_info = bdd_manager_get_info(bddm, t);
  bv_term_compute_value(terms, t, (bvconstant_t**) children_values.data, &t_info->value);

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: new value for ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), t_info->value.data, t_info->value.bitsize);
    ctx_trace_printf(bddm->ctx, "\n");
  }

  // Remove temp
  delete_pvector(&children_values);
}

static inline
void bdd_manager_compute_bdd(bdd_manager_t* bddm, term_t t) {

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: computing BDD for ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "unassigned variable: ");
    ctx_trace_term(bddm->ctx, bddm->unassigned_var);
  }

  uint32_t i;
  term_t t_i;
  BDD** bdds_i;
  term_table_t* terms = bddm->ctx->terms;

  // We store children BDD pointers here
  pvector_t children_bdds;
  init_pvector(&children_bdds, 0);

  // Term data
  term_info_t* t_info = bdd_manager_get_info(bddm, t);

  // Remove all previous BDDs and make a new ones
  bdd_manager_delete_vec(bddm, t_info->v);

  // Negation
  if (is_neg_term(t)) {
    t_i = unsigned_term(t);
    bdds_i = bdd_manager_get_bdds_same_size(bddm, t_i);
    pvector_push(&children_bdds, bdds_i);
  } else {
    term_kind_t t_kind = term_kind(terms, t);

    // First, get all the children BDDs
    switch (t_kind) {
    case OR_TERM:
    case EQ_TERM:
    case BV_EQ_ATOM:
    case BV_GE_ATOM:
    case BV_SGE_ATOM:
    case BV_ARRAY:
    case BV_DIV:
    case BV_REM:
    case BV_SDIV:
    case BV_SREM:
    case BV_SMOD:
    case BV_SHL:
    case BV_LSHR:
    case BV_ASHR: {
      // Get the children BDDs
      composite_term_t* t_comp = composite_term_desc(terms, t);
      for (i = 0; i < t_comp->arity; ++i) {
        t_i = t_comp->arg[i];
        bdds_i = bdd_manager_get_bdds_same_size(bddm, t_i);
        pvector_push(&children_bdds, bdds_i);
      }
      break;
    }
    case BIT_TERM:
      t_i = bit_term_arg(terms, t);
      bdds_i = bdd_manager_get_bdds_same_size(bddm, t_i);
      pvector_push(&children_bdds, bdds_i);
      break;
    case BV_POLY: {
      bvpoly_t* t_poly = bvpoly_term_desc(terms, t);
      for (uint32_t i = 0; i < t_poly->nterms; ++i) {
        if (t_poly->mono[i].var == const_idx) continue;
        t_i = t_poly->mono[i].var;
        bdds_i = bdd_manager_get_bdds_same_size(bddm, t_i);
        pvector_push(&children_bdds, bdds_i);
      }
      break;
    }
    case BV64_POLY: {
      // Polynomial, allocate data for the subterms (variables)
      bvpoly64_t* t_poly = bvpoly64_term_desc(terms, t);
      for (uint32_t i = 0; i < t_poly->nterms; ++i) {
        if (t_poly->mono[i].var == const_idx) continue;
        t_i = t_poly->mono[i].var;
        bdds_i = bdd_manager_get_bdds_same_size(bddm, t_i);
        pvector_push(&children_bdds, bdds_i);
      }
      break;
    }
    case POWER_PRODUCT: {
      pprod_t* t_pprod = pprod_term_desc(terms, t);
      for (uint32_t i = 0; i < t_pprod->len; ++ i) {
        t_i = t_pprod->prod[i].var;
        bdds_i = bdd_manager_get_bdds_same_size(bddm, t_i);
        pvector_push(&children_bdds, bdds_i);
      }
      break;
    }
    default:
      // Shouldn't be here
      assert(false);
    }
  }

  // We have the children values compute
  t_info->v = bdds_compute_bdds(bddm->cudd, terms, t, &children_bdds);

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: BDD done for ");
    ctx_trace_term(bddm->ctx, t);
  }

  // Remove temp
  delete_pvector(&children_bdds);
}

static
bddvec_t bdd_manager_get_term_bdds(bdd_manager_t* bddm, term_t t, uint32_t bitsize) {

  uint32_t i;
  term_t t_recompute;

  assert(bddm->visited.nelems == 0);
  assert(bddm->bdd_recompute.size == 0);
  assert(bddm->value_recompute.size == 0);

  // Make sure we have allocated the data for t
  bdd_manager_ensure_term_data(bddm, t, bitsize);

  // Recompute timestamps
  uint32_t bdd_timestamp = 0, value_timestamp = 0;
  bool contains_unassigned = bdd_manager_recompute_timestamps(bddm, t, bitsize, &bdd_timestamp, &value_timestamp);
  (void) contains_unassigned;
  assert(contains_unassigned);

  // Recompute needed values
  for (i = 0; i < bddm->value_recompute.size; ++ i) {
    t_recompute = bddm->value_recompute.data[i];
    bdd_manager_compute_value(bddm, t_recompute);
  }

  // Recompute needed BDDs
  for (i = 0; i < bddm->bdd_recompute.size; ++ i) {
    assert(i == 0 || t_recompute != bddm->bdd_recompute.data[i]);
    t_recompute = bddm->bdd_recompute.data[i];
    bdd_manager_compute_bdd(bddm, t_recompute);
  }

  // Clear the cache
  int_hmap_reset(&bddm->visited);
  ivector_reset(&bddm->bdd_recompute);
  ivector_reset(&bddm->value_recompute);

  term_info_t* t_info = bdd_manager_get_info(bddm, t);

  // Return the BDDs
  return t_info->v;
}

void bdd_manager_set_bv_value(bdd_manager_t* bddm, term_t t, const bvconstant_t* value) {
  assert(bv_term_get_type(bddm->ctx->terms, t) == BV_TERM_VARIABLE);
  assert(bv_term_bitsize(bddm->ctx->terms, t) == value->bitsize);

  term_info_t* t_info = bdd_manager_get_info(bddm, t);
  assert(t_info->value.bitsize == value->bitsize);

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: setting bit value for ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), value->data, value->bitsize);
    ctx_trace_printf(bddm->ctx, "\nold_value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), t_info->value.data, t_info->value.bitsize);
    ctx_trace_printf(bddm->ctx, "\n");
  }

  if (t_info->value_timestamp == 0 || !bvconst_eq(t_info->value.data, value->data, value->width)) {
    bvconstant_copy(&t_info->value, value->bitsize, value->data);
    t_info->value_timestamp = ++ bddm->timestamp;
  }

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "set value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), t_info->value.data, t_info->value.bitsize);
    ctx_trace_printf(bddm->ctx, "\n");
  }
}

void bdd_manager_set_bool_value(bdd_manager_t* bddm, term_t t, bool value) {
  assert(bv_term_get_type(bddm->ctx->terms, t) == BV_TERM_VARIABLE);
  assert(bv_term_bitsize(bddm->ctx->terms, t) == 1);

  term_info_t* t_info = bdd_manager_get_info(bddm, t);

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: setting bit value for ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "value = %s\n", (value ? "true" : "false"));
    ctx_trace_printf(bddm->ctx, "old_value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), t_info->value.data, t_info->value.bitsize);
    ctx_trace_printf(bddm->ctx, "\n");
  }

  assert(t_info->value.bitsize == 1);

  if (t_info->value_timestamp == 0 || bvconst_tst_bit(t_info->value.data, 0) != value) {
    if (value) {
      bvconst_set_bit(t_info->value.data, 0);
    } else {
      bvconst_clr_bit(t_info->value.data, 0);
    }
    t_info->value_timestamp = ++ bddm->timestamp;
  }

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "set value = ");
    bvconst_print(ctx_trace_out(bddm->ctx), t_info->value.data, t_info->value.bitsize);
    ctx_trace_printf(bddm->ctx, "\n");
  }
}

bddvec_t bdd_manager_get(bdd_manager_t* bddm, term_t t, term_t x) {

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: creating BDD for ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "unit in variable ");
    ctx_trace_term(bddm->ctx, x);
  }

  // Set the unit variable
  assert(bddm->unassigned_var == NULL_TERM);
  bddm->unassigned_var = x;

  // Compute
  assert(is_boolean_term(bddm->ctx->terms, t));
  bddvec_t result = bdd_manager_get_term_bdds(bddm, t, 1);

  // Unset the unit variable
  bddm->unassigned_var = NULL_TERM;

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: BDD for ");
    ctx_trace_term(bddm->ctx, t);
    ctx_trace_printf(bddm->ctx, "unit in variable ");
    ctx_trace_term(bddm->ctx, x);
    bdd_manager_print_bddvec(bddm, result, ctx_trace_out(bddm->ctx));
    ctx_trace_printf(bddm->ctx, "\n");
  }

  return result;
}

BDD** bdd_manager_get_bdds(const bdd_manager_t* bddm, bddvec_t v) {
  return bddvec_manager_get_bdds(&bddm->bdds, v);
}

void bdd_manager_detach(bdd_manager_t* bddm, bddvec_t v) {
  if (v.id == bddvec_null_id) return;
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  bdds_detach(bddm->cudd, v_bdds, v.size);
}

void bdd_manager_delete_vec(bdd_manager_t* bddm, bddvec_t v) {
  if (v.id == bddvec_null_id) return;
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  bdds_clear(bddm->cudd, v_bdds, v.size);
  bddvec_manager_delete_vec(&bddm->bdds, v);
}

bddvec_t bdd_manager_new_copy(bdd_manager_t* bddm, bddvec_t v) {
  if (v.id == bddvec_null_id) return v;
  bddvec_t result = bddvec_manager_new_vec(&bddm->bdds, v.size);
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  BDD** result_bdds = bddvec_manager_get_bdds(&bddm->bdds, result);
  for (uint32_t i = 0; i < v.size; ++ i) {
    result_bdds[i] = v_bdds[i];
  }
  bdds_attach(result_bdds, v.size);
  return result;
}

bddvec_t bdd_manager_new_vec(bdd_manager_t* bddm, uint32_t size) {
  return bddvec_manager_new_vec(&bddm->bdds, size);
}

void bdd_manager_attach(bdd_manager_t* bddm, bddvec_t v) {
  if (v.id == bddvec_null_id) return;
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  bdds_attach(v_bdds, v.size);
}

void bdd_manager_print_bddvec(const bdd_manager_t* bddm, bddvec_t v, FILE* out) {
  BDD** bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  bdds_print(bddm->cudd, bdds, v.size, out);
}

bool bdd_manager_is_empty(const bdd_manager_t* bddm, bddvec_t v) {
  // Check that BDD != false, i.e., there is an assignment so that it is true
  assert(v.id != bddvec_null_id);
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  if (v.size == 1 && v_bdds[0] == bddm->bdd_false) {
    return true;
  }
  return false;
}

bool bdd_manager_is_point(const bdd_manager_t* bddm, bddvec_t v, uint32_t bitsize) {

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, "bdd_manager: checking if point of size %d: ", bitsize);
    bdd_manager_print_bddvec(bddm, v, ctx_trace_out(bddm->ctx));
    ctx_trace_printf(bddm->ctx, "\n");
  }

  // Check if BDD is a point: has to be a cube
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  bool result =  bdds_is_point(bddm->cudd, v_bdds, v.size, bitsize);

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd")) {
    ctx_trace_printf(bddm->ctx, result ? "yes\n" : "no\n");
  }

  return result;
}

/**
 * We're intersecting two vectors of BDDs. Each element of the vector is over
 * disjoint variables in v1 and v2.
 *
 * Inputs are attached BDDs.
 * Return value is always freshly allocated attached.
 */
bddvec_t bdd_manager_intersect(bdd_manager_t* bddm, bddvec_t v1, bddvec_t v2) {

  CUDD* cudd = bddm->cudd;

  BDD** v1_bdds = bddvec_manager_get_bdds(&bddm->bdds, v1);
  BDD** v2_bdds = bddvec_manager_get_bdds(&bddm->bdds, v2);

  pvector_t* result_bdds_tmp = &bddm->bdd_temp;
  assert(result_bdds_tmp->size == 0);

  // Just copy over v1 BDDs (and attach, we're copying)
  bool is_false = false;
  for (uint32_t i = 0; !is_false && i < v1.size; ++ i) {
    BDD* bdd = v1_bdds[i];
    if (bdd == bddm->bdd_false) {
      assert(v1.size == 1); // can be only false, we still add but we're done
      is_false = true;
    }
    assert(bdd != bddm->bdd_true || v1.size == 1); // can be only true, we still add
    pvector_push(result_bdds_tmp, bdd);
  }
  bdds_attach((BDD**) result_bdds_tmp->data, result_bdds_tmp->size);

  if (!is_false) {
    bdds_disjoint_set_add(cudd, v2_bdds, v2.size, result_bdds_tmp);
  }

  // Allocate the BDD vector and copy over the BDDs
  assert(result_bdds_tmp->size > 0);
  bddvec_t result = bddvec_manager_new_vec(&bddm->bdds, result_bdds_tmp->size);
  BDD** result_bdds = bddvec_manager_get_bdds(&bddm->bdds, result);
  for (uint32_t i = 0; i < result_bdds_tmp->size; ++ i) {
    assert(result_bdds_tmp->data[i] != NULL);
    result_bdds[i] = result_bdds_tmp->data[i];
    assert(result_bdds[i] != bddm->bdd_true || result_bdds_tmp->size == 1);
    assert(result_bdds[i] != bddm->bdd_false|| result_bdds_tmp->size == 1);
  }

  if (ctx_trace_enabled(bddm->ctx, "mcsat::bv::bdd::intersect")) {
    ctx_trace_printf(bddm->ctx, "intersect result:\n");
    bdd_manager_print_bddvec(bddm, result, ctx_trace_out(bddm->ctx));
    ctx_trace_printf(bddm->ctx, "\n");
  }

  pvector_reset(result_bdds_tmp);

  return result;
}

void bdd_manager_pick_value(bdd_manager_t* bddm, term_t x, bddvec_t v, bvconstant_t* out) {
  term_info_t* x_info = bdd_manager_get_info(bddm, x);
  const bvconstant_t* prev_value = &x_info->value;
  BDD** x_bdds = bddvec_manager_get_bdds(&bddm->bdds, x_info->v);
  assert(out->bitsize == prev_value->bitsize);
  if (bdd_manager_is_model(bddm, x, v, prev_value)) {
    // If cached value works, just use it
    bvconstant_copy(out, out->bitsize, prev_value->data);
  } else {
    // Otherwise, get individual model segments
    uint32_t i;
    BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
    bvconstant_set_all_zero(out, out->bitsize);
    for (i = 0; i < v.size; ++ i) {
      bdds_get_model(bddm->cudd, x_bdds, v_bdds[i], out);
    }
  }
}

bool bdd_manager_is_model(bdd_manager_t* bddm, term_t x, bddvec_t v, const bvconstant_t* x_value) {
  uint32_t i;
  term_info_t* x_info = bdd_manager_get_info(bddm, x);
  BDD** x_bdds = bddvec_manager_get_bdds(&bddm->bdds, x_info->v);
  BDD** v_bdds = bddvec_manager_get_bdds(&bddm->bdds, v);
  for (i = 0; i < v.size; ++ i) {
    if (!bdds_is_model(bddm->cudd, x_bdds, v_bdds[i], x_value))
      return false;
  }
  return true;
}

