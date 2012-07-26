/*
 * COMPUTE THE VALUE OF A TERM IN A MODEL
 */

#ifndef __MODEL_EVAL_H
#define __MODEL_EVAL_H

#include <stdint.h>
#include <setjmp.h>

#include "models.h"
#include "int_stack.h"
#include "int_hash_map.h"


/*
 * Error codes returned by eval_in_model
 * - all are negative integers
 * - a non-negative value means no error
 * - the null_value is -1 so we start from -2
 */
enum {
  MDL_EVAL_INTERNAL_ERROR = -2,
  MDL_EVAL_UNKNOWN_TERM = -3,
  MDL_EVAL_FREEVAR_IN_TERM = -4,
  MDL_EVAL_QUANTIFIER = -5,
  MDL_EVAL_FAILED = -6, // function equality involved
};



/*
 * Evaluator structure:
 * - pointer to a model + term_table + value_table 
 *   (term table and value table are extracted from
 *    model when the evaluator is initialized)
 * - cache: keeps track of the value of evaluated terms
 * - env: jump buffer for error handling
 * - stack of integer arrays
 */
typedef struct evaluator_s {
  model_t *model;
  term_table_t *terms;
  value_table_t *vtbl;
  int_hmap_t cache;
  int_stack_t stack;
  jmp_buf env;
} evaluator_t;




/*
 * Initialization for the given model
 */
extern void init_evaluator(evaluator_t *eval, model_t *model);


/*
 * Deletion: free all memory
 */
extern void delete_evaluator(evaluator_t *eval);


/*
 * Reset: empty the cache and delete all temporary objects
 * created in eval->model.vtbl
 */
extern void reset_evaluator(evaluator_t *eval);


/*
 * Compute the value of term t in the model
 * - t must be a valid term
 * - return a negative code if there's an error
 * - return the id of a concrete objects of eval->model.vtbl
 *
 * Evaluation may create new objects. All these new objects are
 * marked as temporary objects and can be deleted by calling
 * reset_evaluator.
 */
extern value_t eval_in_model(evaluator_t *eval, term_t t);



#endif /* __MODEL_EVAL_H */
