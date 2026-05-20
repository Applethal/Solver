#include "core.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
// TODO: Work on the objective function, btw the objective function already has the constant from variables that hit the upper bound so just take the RHS values 

#define MIN3IDX(a,b,c) (((a)<=(b)&&(a)<=(c)) ? 1 : (((b)<(a)&&(b)<=(c)) ? 2 : 3))


void Get_ObjectiveFunctionBounded(Model *model, double *rhs_vector) {
  size_t n = model->num_constraints;

  }


void BoundedSimplex(Model *model) {

  printf("=== MODEL INPUT DEBUG ===\n");
printf("RHS: ");
for (int i = 0; i < model->num_constraints; i++)
    printf("%f ", model->rhs_vector[i]);
printf("\n");

for (int i = 0; i < model->num_vars; i++)
    printf("x%d: ub=%f originallb=%f value=%f status=%d\n",
           i, model->coeffs[i].ub, model->coeffs[i].originallb,
           model->coeffs[i].value, model->coeffs[i].status);

printf("LHS matrix:\n");
int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;
for (int i = 0; i < model->num_constraints; i++) {
    for (int j = 0; j < total_cols; j++)
        printf("%6.2f ", model->lhs_matrix[i][j]);
    printf("\n");
}
printf("Objective constant: %lf\n", model->objective_function);

printf("=========================\n");

  int termination = 0;
  size_t n = model->num_constraints;
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) * 10 + 1;

  while (termination != 1) {

    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      return;
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
      double reduced_cost =
          Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
      if (reduced_cost > best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0) {
        feasibility_check++;
      }
    }
    printf("Best reduced cost is %lf, entering variable is x%i which has the index %i \n", best_reduced_cost, entering_var+1, entering_var_idx);
    if (feasibility_check == model->non_basics_count) {
      printf("Optimal solution found! \n");
      //Get_ObjectiveFunctionBounded(model, original_RHS);
      termination++;
      for (size_t i = 0; i < n; i++) {
        free(B_inv[i]);
      }
      free(B_inv);
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
      if (Pivot[i] <= 1e-6) {
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
      for (size_t i = 0; i < n; i++) {
        free(B_inv[i]);
      }
      free(B_inv);
      free(Simplex_multiplier);

      return;
    }
  int winner_ratio = MIN3IDX(best_ratio, best_ratio_neg, model->coeffs[entering_var].ub); // 1 for ratio 1, 2 for ratio 2 and 3 if upper bound wins 
  switch (winner_ratio) {
  case 1:
    printf("Ratio 1 wins here, exiting var is x%i and has an index of %i\n", exiting_var+1, exiting_var_idx);
    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;

    break;
  case 2:
    // code block
    printf("Ratio 2 wins here, exiting var is x%i and has an index of %i\n", exiting_var_neg+1, exiting_var_idx_neg);
    for (size_t i = 0; i < n; i++) {
        model->rhs_vector[i] -= model->lhs_matrix[i][exiting_var_neg] * model->coeffs[exiting_var_neg].ub;
        model->lhs_matrix[i][exiting_var_neg] *= -1;
    }
    model->coeffs[exiting_var_neg].value *= -1;
    model->basics_vector[exiting_var_idx_neg] = entering_var;
    model->non_basics[entering_var_idx] = exiting_var_neg;

    break;
  case 3:
  printf("Ratio 3 wins here, no exiting variable here, we just substitute the entering variable x%i with its upper bound and flip the entering var's sign\n", entering_var +1);

  // Not sure if this is the right thing to do:
  for (size_t i = 0; i < n ; i++) {
      model->rhs_vector[i] -= model->lhs_matrix[i][entering_var] * model->coeffs[entering_var].ub;
      model->lhs_matrix[i][entering_var] *= -1; 
        
  }
  model->objective_function += model->coeffs[entering_var].ub * model->coeffs[entering_var].value; 
  model->coeffs[entering_var].value *= -1;
  
  
}
    
  printf("=== MODEL INPUT DEBUG ===\n");
printf("RHS: ");
for (int i = 0; i < model->num_constraints; i++)
    printf("%f ", model->rhs_vector[i]);
printf("\n");

for (int i = 0; i < model->num_vars; i++)
    printf("x%d: ub=%f originallb=%f value=%f status=%d\n",
           i, model->coeffs[i].ub, model->coeffs[i].originallb,
           model->coeffs[i].value, model->coeffs[i].status);

printf("LHS matrix:\n");
int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;
for (int i = 0; i < model->num_constraints; i++) {
    for (int j = 0; j < total_cols; j++)
        printf("%6.2f ", model->lhs_matrix[i][j]);
    printf("\n");
}
printf("Objective constant: %lf\n", model->objective_function);

printf("Basics vector:\n");
for (int i = 0; i < model->non_basics_count; i++) {
  printf(" %i ", model->non_basics[i]+1);

}
printf("\n non basics: \n");

for (int i = 0; i < 2; i++) {

  printf(" %i ", model->basics_vector[i]+1 );
}

printf("\n");
printf("=========================\n");


  



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
  TransformBoundedModel(model);
  ValidateModelPointers(model);
  BoundedSimplex(model);
}

void TransformBoundedModel(Model *model) {

  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
  }

  // subtracking the bounds from the RHS values
// TODO: What if RHS becomes negative?  

   for (size_t j = 0; j < model->num_constraints; j++) {
    for (size_t i = 0; i < model->num_vars; i++) {
      model->rhs_vector[j] -= model->lhs_matrix[j][i] * model->coeffs[i].lb;
    }
    if (model->rhs_vector[j] < 0 ) {
        FlipConstraint(model, j); // Yeah this should do, for now. Later I will have to probably destroy the constraint since it becomes redundant if the LHS is positive
    
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
      model->coeffs[col].ub = model->rhs_vector[i];
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
      model->coeffs[surplus_col].ub = model->rhs_vector[i];

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
      model->coeffs[col].status = 1;
    }
  }

  model->non_basics_count = non_basic_idx;
}
