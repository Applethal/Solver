#include "core.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

void Get_ObjectiveFunctionBounded(Model *model, double *rhs_vector) {
    size_t n = model->num_constraints;
    
    for (int i = 0; i < n; i++) {
        int basic_idx = model->basics_vector[i];
        if (model->coeffs[basic_idx].type == ARTIFICIAL && rhs_vector[i] > 1e-6) {
            printf("Model Infeasible. Artificial basic variable has positive value. Terminating!\n");
            return;
        }
    }
    
    // Recover basic variable values
    for (int i = 0; i < n; i++) {
        int basic_idx = model->basics_vector[i];
        if (model->coeffs[basic_idx].type == STANDARD) {
            double original_value = rhs_vector[i] + model->coeffs[basic_idx].originallb;

            model->objective_function += model->coeffs[basic_idx].value * rhs_vector[i];

            printf("Value of variable x%i is %f\n", basic_idx, original_value);
        }
    }
    
    // Recover non-basic variables at upper bound
    for (int i = 0; i < model->non_basics_count; i++) {
        int non_basic_idx = model->non_basics[i];
        if (model->coeffs[non_basic_idx].status == UPPER && model->coeffs[non_basic_idx].type == STANDARD) {

            double original_value = model->coeffs[non_basic_idx].ub +  model->coeffs[non_basic_idx].originallb;

            model->objective_function += model->coeffs[non_basic_idx].value *  model->coeffs[non_basic_idx].ub;
            printf("Value of variable x%i is %f\n", non_basic_idx, original_value);
        }
    }
    
    model->objective_function += model->ObjectiveConstant;
    
    if (model->objective == 0) {
        model->objective_function *= -1;
    }
    
}
 
 
void BoundedSimplex(Model *model) {
    
    // For testing only, I don't have a parser yet

    model->rhs_vector[0] = 14;
model->rhs_vector[1] = 7;
model->rhs_vector[2] = 14;
model->rhs_vector[3] = 11;
model->rhs_vector[4] = 9;

model->coeffs[0].status = LOWER;
model->coeffs[1].status = LOWER;
model->coeffs[2].status = LOWER;
model->coeffs[3].status = LOWER;
model->coeffs[4].status = LOWER;

model->coeffs[0].ub = 5;
model->coeffs[1].ub = 5;
model->coeffs[2].ub = 4;
model->coeffs[3].ub = 5;
model->coeffs[4].ub = 3;

model->ObjectiveConstant = 13;

model->coeffs[0].originallb = 0;
model->coeffs[1].originallb = 1;
model->coeffs[2].originallb = 0;
model->coeffs[3].originallb = 2;
model->coeffs[4].originallb = 0;
 
    // Big-M penalty applied AFTER negation
    int artificial_start = model->num_vars + model->slacks_surplus_count;
    for (int i = 0; i < model->artificials_count; i++) {
        model->coeffs[artificial_start + i].value =
            -(model->bigM * 2 * model->num_vars);
    }
 
    // Set initial status for all non-basics
    for (int i = 0; i < model->non_basics_count; i++) {
        model->coeffs[model->non_basics[i]].status = LOWER;
    }
 
    int termination = 0;
    size_t n = model->num_constraints;
    int MAX_ITERATIONS = (model->num_vars * model->num_constraints) * 10 + 1;
 
    while (termination != 1) {
 
        
        if (model->solver_iterations == MAX_ITERATIONS) {
            printf("Max iterations reached. Terminating!\n");
            FreeModel(model);
            exit(1);
        }
 
        int feasibility_check = 0;
        printf("Beginning solver iteration %i ... \n", model->solver_iterations);


        double **B_inv = Get_BasisInverse(model, model->solver_iterations);
        double original_RHS[n];
        for (size_t i = 0; i < n; i++) {
            original_RHS[i] = model->rhs_vector[i];
        }
 
        double *Simplex_multiplier = Get_SimplexMultiplier(model, B_inv);
 
        int entering_var_idx = -1;
        int entering_var = -1;
        UpdateRhs(model, original_RHS, B_inv);
 
        double best_reduced_cost = 0;
        for (size_t i = 0; i < model->non_basics_count; i++) {
            int non_basic_idx = model->non_basics[i];
            double reduced_cost = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
 
            if (model->coeffs[non_basic_idx].status == LOWER && reduced_cost <= 0) {
                feasibility_check++;
            } else if (model->coeffs[non_basic_idx].status == UPPER && reduced_cost >= 0) {
                feasibility_check++;
            }
 
            if (model->coeffs[non_basic_idx].status == LOWER && reduced_cost > best_reduced_cost) {
                best_reduced_cost = reduced_cost;
                entering_var_idx = i;
                entering_var = non_basic_idx;
            } else if (model->coeffs[non_basic_idx].status == UPPER && reduced_cost < 0 && fabs(reduced_cost) > best_reduced_cost) {
                best_reduced_cost = fabs(reduced_cost);
                entering_var_idx = i;
                entering_var = non_basic_idx;
            }
        }
 
 
        if (feasibility_check == model->non_basics_count) {
            printf("Optimal solution found! \n");
            Get_ObjectiveFunctionBounded(model, original_RHS);
            termination++;
            for (size_t i = 0; i < n; i++) {
                free(B_inv[i]);
            }
            free(B_inv);
            free(Simplex_multiplier);
            break;
        }
 
        double *Pivot = Get_pivot_column(B_inv, model, entering_var);
 
        // Determine direction: UPPER-status entering var moves downward
        double direction = (model->coeffs[entering_var].status == UPPER) ? -1.0 : 1.0;
 
        int exiting_var_idx = -1;
        int exiting_var = -1;
        double best_ratio = DBL_MAX;
        double ratio;
 
        for (size_t i = 0; i < n; i++) {
            ratio = DBL_MAX;
            double piv = Pivot[i] * direction;
 
            if (piv > 1e-6) {
                // Basic variable hits its lower bound (0)
                ratio = original_RHS[i] / piv;
            } else if (piv < -1e-6 && model->coeffs[model->basics_vector[i]].ub < DBL_MAX) {
                // Basic variable hits its upper bound
                ratio = (model->coeffs[model->basics_vector[i]].ub - original_RHS[i]) / -piv;
            }
 
            if (ratio < best_ratio) {
                best_ratio = ratio;
                exiting_var_idx = i;
                exiting_var = model->basics_vector[i];
            }
        }
 
        if (model->coeffs[entering_var].ub <= best_ratio) {
            // Entering variable flips to/from its bound without entering the basis
 
            best_ratio = model->coeffs[entering_var].ub;
 
            if (model->coeffs[entering_var].status == LOWER) {
                model->coeffs[entering_var].status = UPPER;
                // Subtract the column contribution from stored RHS
                for (size_t i = 0; i < n; i++) {
                    model->rhs_vector[i] -= model->lhs_matrix[i][entering_var]
                                            * model->coeffs[entering_var].ub;
                }
            } else {
                // Was UPPER, flip back to LOWER
                model->coeffs[entering_var].status = LOWER;
                for (size_t i = 0; i < n; i++) {
                    model->rhs_vector[i] += model->lhs_matrix[i][entering_var]
                                            * model->coeffs[entering_var].ub;
                }
            }
 
        } else {
            model->non_basics[entering_var_idx] = exiting_var;
            model->basics_vector[exiting_var_idx] = entering_var;
            model->coeffs[entering_var].status = BASIC;
 
            // Exiting variable goes to lower or upper bound depending on which it hit
            double piv = Pivot[exiting_var_idx] * direction;
            if (piv > 1e-6) {
                // Hit lower bound
                model->coeffs[exiting_var].status = LOWER;
            } else {
                // Hit upper bound — update stored RHS for this new UPPER non-basic
                model->coeffs[exiting_var].status = UPPER;
                for (size_t i = 0; i < n; i++) {
                    model->rhs_vector[i] -= model->lhs_matrix[i][exiting_var]
                                            * model->coeffs[exiting_var].ub;
                }
            }
        }
 
 
        model->solver_iterations++;
        free(Pivot);
        for (size_t i = 0; i < n; i++) {
            free(B_inv[i]);
        }
        free(B_inv);
        free(Simplex_multiplier);
    }
}

void SolveBounded(Model *model) {
  // I need to remember that I won't adjust bounds here at all, bound
  // arithmetics will have to be considered while working on the algorithm The
  // reason for this is because I want to save some CPU cycles
  TransformModel(model);
  ValidateModelPointers(model);
  BoundedSimplex(model);
  
}

void TransformBoundedModel(Model *model) {

  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
  }

  // subtracking the bounds from the RHS values

  for (size_t j = 0; j < model->num_constraints; j++) {
    for (size_t i = 0; i < model->num_vars; i++) {
      model->rhs_vector[j] -= model->lhs_matrix[j][i] * model->coeffs[i].lb;
    }
  }

  // Pushing LBs to 0, reduce the range of the UBs and add a constant to
  // objective function to consider these changes

  for (size_t i = 0; i < model->num_vars; i++) {

    model->ObjectiveConstant += model->coeffs[i].value * model->coeffs[i].lb;
    model->coeffs[i].ub -= model->coeffs[i].lb;
    model->coeffs[i].originallb = model->coeffs[i].lb;
    model->coeffs[i].lb = 0;
  }

  int total_cols =
      model->num_vars + model->slacks_surplus_count + model->artificials_count;

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
      model->coeffs[col].constraint_idx = 0;
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
      model->coeffs[surplus_col].constraint_idx = 0;
      model->coeffs[surplus_col].lb = 0;
      model->coeffs[surplus_col].ub = DBL_MAX;

      model->coeffs[artif_col].type = ARTIFICIAL;
      model->coeffs[artif_col].value = -1.0;
      model->coeffs[artif_col].constraint_idx = 0;
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
      model->coeffs[artif_col].value = -1.0;
      model->coeffs[artif_col].constraint_idx = 0;
      model->coeffs[artif_col].lb = 0;
      model->coeffs[artif_col].ub = DBL_MAX;

      model->artificials_vector[artificial_idx] = artif_col;
      model->basics_vector[i] = artif_col;

      artificial_idx++;
    }
  }

  // Artificial values will get a value of 1. In previous versions, I gave them
  // a value of max(coeffs) * 2 to penalize the objective function, now I can
  // track them with their types Convert MINIMIZE to MAXIMIZE by negating all
  // coefficients
  if (model->objective == 0) {
    for (int i = 0; i < total_cols; i++)
      model->coeffs[i].value *= -1;
  }

  // Build non-basics list
  model->non_basics = (int *)malloc(model->num_vars * sizeof(int));
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
