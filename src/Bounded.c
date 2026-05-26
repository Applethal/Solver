#include "core.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define REINVERSION_FREQ 50
#define MIN3IDX(a,b,c) (((a)<=(b)&&(a)<=(c)) ? 1 : (((b)<(a)&&(b)<=(c)) ? 2 : 3))

// printing part was vibecoded
void Get_ObjectiveFunctionBounded(Model *model, double *original_RHS) {
    size_t n = model->num_constraints;
 
    printf("\n%-12s %-10s %-16s %-20s\n", "Variable", "Value", "Value bounds", "Status");
    printf("%-12s %-10s %-16s %-20s\n",   "--------", "-----", "------------", "------");
 
    // basic variables
    for (int i = 0; i < (int)n; i++) {
        int basic_var = model->basics_vector[i];
        if (model->coeffs[basic_var].type != STANDARD) continue;
 
        double c_j = (model->coeffs[basic_var].status == UPPER)
                     ? model->coeffs[basic_var].value * -1
                     : model->coeffs[basic_var].value;
        double x_j = (model->coeffs[basic_var].status == UPPER)
                     ? model->coeffs[basic_var].ub - original_RHS[i]
                     : original_RHS[i];
 
        model->objective_function += (c_j * x_j) * model->objective;
 
        double lb = model->coeffs[basic_var].originallb;
        double ub = model->coeffs[basic_var].ub + lb;
        x_j += lb;
 
        char bounds_str[32];
        snprintf(bounds_str, sizeof(bounds_str), "[%.4g, %.4g]", lb, ub);
        printf("x%-11i %-10.4f %-16s %-20s\n", basic_var+1, x_j, bounds_str, "Basic");
    }
 
    // non-basic variables
    for (int i = 0; i < model->non_basics_count; i++) {
        int non_basic = model->non_basics[i];
        if (model->coeffs[non_basic].type != STANDARD) continue;
 
        double lb = model->coeffs[non_basic].originallb;
        double ub = model->coeffs[non_basic].ub + lb;
 
        char bounds_str[32];
        snprintf(bounds_str, sizeof(bounds_str), "[%.4g, %.4g]", lb, ub);
 
        if (model->coeffs[non_basic].status == UPPER) {
            printf("x%-11i %-10.4f %-16s %-20s\n", non_basic+1, ub, bounds_str, "At upper bound");
        } else {
            printf("x%-11i %-10.4f %-16s %-20s\n", non_basic+1, lb, bounds_str, "At lower bound");
        }
    }
 
    printf("\nObjective function: %f\n", model->objective_function);
}
 


void BoundedSimplex(Model *model) {



  int termination = 0;
  size_t n = model->num_constraints;
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) * 10 + 1;
  double **B_inv = Get_BasisInverse(model, model->solver_iterations);

  while (termination != 1) {

    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      break;
    }

    printf("Beginning solver iteration %i ... \n", model->solver_iterations);

    
    double original_RHS[n];
    for (size_t i = 0; i < n; i++) {
      original_RHS[i] = model->rhs_vector[i];
    }

    double *Simplex_multiplier = Get_SimplexMultiplier(model, B_inv);
    

  
    int entering_var_idx = -1;
    int entering_var = -1;
    UpdateRhs(model, original_RHS, B_inv);
 

    double best_reduced_cost = -DBL_MAX;
    for (size_t i = 0; i < model->non_basics_count; i++) {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
      if (reduced_cost > best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
    }
    if (best_reduced_cost <= 1e-9) {
      printf("Optimal solution found! \n");
      Get_ObjectiveFunctionBounded(model, original_RHS);
      termination++;
      
      free(Simplex_multiplier);
      break;
    }

    double *Pivot = Get_pivot_column(B_inv, model, entering_var);

    
    int exiting_var_idx = -1;
    int exiting_var = -1; 
    double best_ratio = DBL_MAX;
    double ratio;
    // neg here refers to a_{ij} that is negative, if there is none, then this ratio test is best left at infinity
    double best_ratio_neg = DBL_MAX;
    double ratio_neg = DBL_MAX; 
    int exiting_var_idx_neg = -1;
    int exiting_var_neg = -1;

    for (size_t i = 0; i < n; i++) {
      if (Pivot[i] < -1e-6) {
          int idx = model->basics_vector[i]; 
          ratio_neg = (model->coeffs[idx].ub - original_RHS[i]) / -Pivot[i];
          
          if (ratio_neg < best_ratio_neg) {
              best_ratio_neg = ratio_neg;
              exiting_var_idx_neg = i; 
              exiting_var_neg = model->basics_vector[i];
          }
}
      if (Pivot[i] >= 1e-6) {
          
          ratio = original_RHS[i] / Pivot[i];
          if (ratio < best_ratio) {
              best_ratio = ratio;
              exiting_var_idx = i;
              exiting_var = model->basics_vector[i];    
          }
      }
}
    if (exiting_var_idx == -1 && exiting_var_idx_neg == -1) { 
      printf("LP is unbounded! Terminating!\n");
      free(Pivot);
      
      free(Simplex_multiplier);

      break;
    }
  int winner_ratio = MIN3IDX(best_ratio, best_ratio_neg, model->coeffs[entering_var].ub); // 1 for ratio 1, 2 for ratio 2 and 3 if upper bound wins, bless macros! 
  switch (winner_ratio) {
  case 1:
    Update_BasisInverse(B_inv, Pivot, exiting_var_idx, n);

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;
    model->coeffs[exiting_var].status = LOWER;
    model->coeffs[entering_var].status = LOWER;

    break;
  case 2:
    Update_BasisInverse(B_inv, Pivot, exiting_var_idx_neg, n);
    for (size_t i = 0; i < n; i++) {
        model->rhs_vector[i] -= model->lhs_matrix[i][exiting_var_neg] * model->coeffs[exiting_var_neg].ub;
        model->lhs_matrix[i][exiting_var_neg] *= -1;
    }
    model->objective_function += (model->coeffs[exiting_var_neg].ub * model->coeffs[exiting_var_neg].value) * model->objective; 
    model->coeffs[exiting_var_neg].value *= -1;

    model->basics_vector[exiting_var_idx_neg] = entering_var;
    model->non_basics[entering_var_idx] = exiting_var_neg;
    model->coeffs[exiting_var_neg].status = UPPER;
    model->coeffs[entering_var].status = LOWER;

    break;
  case 3:

  for (size_t i = 0; i < n ; i++) {
      model->rhs_vector[i] -= model->lhs_matrix[i][entering_var] * model->coeffs[entering_var].ub;
      model->lhs_matrix[i][entering_var] *= -1; 
        
  }
  model->coeffs[entering_var].status = UPPER;
  model->objective_function += (model->coeffs[entering_var].ub * model->coeffs[entering_var].value) * model->objective;   
  model->coeffs[entering_var].value *= -1;
  
  
}
 
if (model->solver_iterations % REINVERSION_FREQ == 0) {
    for (int i = 0; i < n; i++) free(B_inv[i]);
    free(B_inv);
    B_inv = Get_BasisInverse(model, model->solver_iterations);
}


    model->solver_iterations++;

    
    free(Pivot);
        free(Simplex_multiplier);
  }
for (size_t i = 0; i < n; i++) {
      free(B_inv[i]);
    }
    free(B_inv);

}

void SolveBounded(Model *model) {
  TransformBoundedModel(model);
  ValidateModelPointers(model);
  BoundedSimplex(model);
}

void TransformBoundedModel(Model *model) {

  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
  }

// TODO: What if RHS becomes negative? For now it is ignored  

   for (size_t j = 0; j < model->num_constraints; j++) {
    for (size_t i = 0; i < model->num_vars; i++) {
      model->rhs_vector[j] -= model->lhs_matrix[j][i] * model->coeffs[i].lb;
    }
  }
  // Pushing LBs to 0, reduce the range of the UBs and add a constant to
  // objective function to consider these changes

  for (size_t i = 0; i < model->num_vars; i++) {

    model->objective_function += (model->coeffs[i].value * model->coeffs[i].lb) * model->objective;
    model->coeffs[i].ub -= model->coeffs[i].lb;
    model->coeffs[i].originallb = model->coeffs[i].lb;
    model->coeffs[i].lb = 0;
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
      // model->coeffs[col].ub = model->rhs_vector[i];
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
      // model->coeffs[surplus_col].ub = model->rhs_vector[i];
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

  model->non_basics = (int *)malloc( (total_cols - model->num_constraints) * sizeof(int));
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
      model->coeffs[col].status = LOWER;
    }
  }

  model->non_basics_count = non_basic_idx;
}
