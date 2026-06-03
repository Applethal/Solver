#include "core.h"
#include "Variable.h"
#include "Constraint.h"


void TransformModel(Model *model) {

  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
  }

  int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;

  size_t memory_needed = (size_t)model->num_constraints * total_cols * sizeof(double);
  size_t max_memory = (size_t)MAX_MEM_BYTES * 1024 * 1024;

  if (memory_needed > max_memory) {
    fprintf(stderr, "Error: Model too large. Requires %zu MB, maximum is %d MB\n",
            memory_needed / (1024 * 1024), MAX_MEM_BYTES);
    exit(1);
  }

  Variable *old_coeffs = model->coeffs;

  // Expand each constraint's lhs_vector to total_cols in-place
  for (int i = 0; i < model->num_constraints; i++) {
    double *old_lhs = model->constraints[i].lhs_vector;
    model->constraints[i].lhs_vector = (double *)calloc(total_cols, sizeof(double));
    for (int j = 0; j < model->num_vars; j++) {
      model->constraints[i].lhs_vector[j] = old_lhs[j];
    }
    free(old_lhs);
  }

  // Expand coefficients array
  model->coeffs = (Variable *)calloc(total_cols, sizeof(Variable));
  for (int i = 0; i < model->num_vars; i++) {
    model->coeffs[i] = old_coeffs[i];
  }
  free(old_coeffs);

  model->artificials_vector = (int *)malloc(model->artificials_count * sizeof(int));
  model->basics_vector      = (int *)calloc(model->num_constraints, sizeof(int));

  int slack_col       = model->num_vars;
  int artificial_col  = model->num_vars + model->slacks_surplus_count;
  int slack_idx       = 0;
  int artificial_idx  = 0;

  for (int i = 0; i < model->num_constraints; i++) {
    if (model->constraints[i].constraints_symbol == 'L') {
      int col = slack_col + slack_idx;
      model->constraints[i].lhs_vector[col] = 1.0;
      model->basics_vector[i]    = col;
      model->coeffs[col].type    = SLACK;
      model->coeffs[col].value   = 0.0;
      model->coeffs[col].lb      = 0;
      model->coeffs[col].ub      = DBL_MAX;
      slack_idx++;

    } else if (model->constraints[i].constraints_symbol == 'G') {
      int surplus_col = slack_col + slack_idx;
      int artif_col   = artificial_col + artificial_idx;

      model->constraints[i].lhs_vector[surplus_col] = -1.0;
      model->constraints[i].lhs_vector[artif_col]   =  1.0;

      model->coeffs[surplus_col].type  = SLACK;
      model->coeffs[surplus_col].value = 0.0;
      model->coeffs[surplus_col].lb    = 0;
      model->coeffs[surplus_col].ub    = DBL_MAX;

      model->coeffs[artif_col].type  = ARTIFICIAL;
      model->coeffs[artif_col].value = -(model->bigM * 2 * model->num_vars);
      model->coeffs[artif_col].lb    = 0;
      model->coeffs[artif_col].ub    = DBL_MAX;

      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      slack_idx++;
      artificial_idx++;

    } else if (model->constraints[i].constraints_symbol == 'E') {
      int artif_col = artificial_col + artificial_idx;
      model->constraints[i].lhs_vector[artif_col] = 1.0;

      model->coeffs[artif_col].type  = ARTIFICIAL;
      model->coeffs[artif_col].value = -(model->bigM * 2 * model->num_vars);
      model->coeffs[artif_col].lb    = 0;
      model->coeffs[artif_col].ub    = DBL_MAX;

      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      artificial_idx++;
    }
  }

  model->non_basics = (int *)malloc((total_cols - model->num_constraints) * sizeof(int));
  int non_basic_idx = 0;
  for (int col = 0; col < total_cols; col++) {
    bool is_basic = false;
    for (int j = 0; j < model->num_constraints; j++) {
      if (model->basics_vector[j] == col) {
        is_basic = true;
        break;
      }
    }
    if (!is_basic) {
      model->non_basics[non_basic_idx++] = col;
    }
  }
  model->non_basics_count = non_basic_idx;
}
