#include "core.h"


void TransformModel(Model *model) {
  // TODO: Make a seperate function for two phase to make it more efficient by directly attributing bigM to the artificials for the BigM method and -1 for the two phase method, maybe as a flag or maybe in the input file
// TODO2: Add bounds to slacks
  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
  }
  int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;

  size_t memory_needed =
      (size_t)model->num_constraints * total_cols * sizeof(double);
  size_t max_memory = (size_t)MAX_MEM_BYTES * 1024 * 1024;

  if (memory_needed > max_memory) {
    fprintf(stderr,
            "Error: Model too large. Requires %zu MB, maximum is %d MB\n",
            memory_needed / (1024 * 1024), MAX_MEM_BYTES);
    exit(1);
  }

  double **old_lhs_matrix = model->lhs_matrix;
  Variable *old_coeffs = model->coeffs;

  // Allocate new expanded lhs_matrix matrix
  model->lhs_matrix =
      (double **)malloc(model->num_constraints * sizeof(double *));
  for (int i = 0; i < model->num_constraints; i++) {
    model->lhs_matrix[i] = (double *)calloc(total_cols, sizeof(double));

    for (int j = 0; j < model->num_vars; j++) {
      model->lhs_matrix[i][j] = old_lhs_matrix[i][j];
    }
  }

  // Allocate new expanded coefficients array
  model->coeffs = (Variable *)calloc(total_cols, sizeof(Variable));

  // Copy original objective coefficients
  for (int i = 0; i < model->num_vars; i++) {
    model->coeffs[i] = old_coeffs[i];
  }

  // Free old arrays
  for (int i = 0; i < model->num_constraints; i++) {
    free(old_lhs_matrix[i]);
  }
  free(old_lhs_matrix);
  free(old_coeffs);

  // Allocate tracking arrays
  model->artificials_vector =
      (int *)malloc(model->artificials_count * sizeof(int));
  model->basics_vector = (int *)calloc(model->num_constraints, sizeof(int));

  // Add slack, surplus, and artificial variables
  int slack_col = model->num_vars;
  int artificial_col = model->num_vars + model->slacks_surplus_count;
  int slack_idx = 0;
  int artificial_idx = 0;

  for (int i = 0; i < model->num_constraints; i++) {
    if (model->constraints_symbols[i] == 'L') {
      // Less-than: add slack variable
      int col = slack_col + slack_idx;
      model->lhs_matrix[i][col] = 1.0;
      model->basics_vector[i] = col;
      model->coeffs[col].type = SLACK;
      model->coeffs[col].value = 0.0;
      model->coeffs[col].lb = 0;
      model->coeffs[col].ub = DBL_MAX;
      slack_idx++;
    } else if (model->constraints_symbols[i] == 'G') {
      // Greater-than: add surplus and artificial
      int surplus_col = slack_col + slack_idx;
      int artif_col = artificial_col + artificial_idx;

      model->lhs_matrix[i][surplus_col] = -1.0;
      model->lhs_matrix[i][artif_col] = 1.0;

      model->coeffs[surplus_col].type = SLACK;
      model->coeffs[surplus_col].value = 0.0;
      model->coeffs[surplus_col].lb = 0;
      model->coeffs[surplus_col].ub = DBL_MAX;

      model->coeffs[artif_col].type = ARTIFICIAL;
      model->coeffs[artif_col].value = -(model->bigM * 2 * model->num_vars);
      model->coeffs[artif_col].lb = 0;
      model->coeffs[artif_col].ub = DBL_MAX;

      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      slack_idx++;
      artificial_idx++;
    } else if (model->constraints_symbols[i] == 'E') {
      // Equality: add artificial variable
      int artif_col = artificial_col + artificial_idx;
      model->lhs_matrix[i][artif_col] = 1.0;

      model->coeffs[artif_col].type = ARTIFICIAL;
      model->coeffs[artif_col].value = -(model->bigM * 2 * model->num_vars);
      model->coeffs[artif_col].lb = 0;
      model->coeffs[artif_col].ub = DBL_MAX;

      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      artificial_idx++;
    }
  }

  //Convert MINIMIZE to MAXIMIZE by negating all coefficients
  // if (model->objective == 0) {
  //   for (int i = 0; i < model->num_vars; i++)
  //     model->coeffs[i].value *= -1;
  // }

  // Build non-basics list
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
