#include "Variable.h"
#include "core.h"
#include <stdio.h>
#include "Constraint.h"


void Solve(Model *model) {
  printf("Warning, Debug mode is disabled. I have not implemented it here yet.");
  TransformModel(model);
  ValidateModelPointers(model);
  if (RunPhase1(model))
    RunPhase2(model);
}

void RunPhase2(Model *model) {

  double solution[model->num_constraints];
  int result = SimplexLoop(model, solution);
  if (result == -1) {
    printf("LP is unbounded!\n");
    return;
  }
  Get_ObjectiveFunction(model, solution);
}

int RunPhase1(Model *model) {
  int total_cols =
    model->num_vars + model->slacks_surplus_count + model->artificials_count;

  // Save real objective coefficients
  double *saved = malloc(total_cols * sizeof(double));
  for (int i = 0; i < total_cols; i++)
    saved[i] = model->coeffs[i].value;

  // Phase 1 objective: -1 for artificials, 0 for everything else
  for (int i = 0; i < total_cols; i++) {
    if (model->coeffs[i].type == ARTIFICIAL)
      model->coeffs[i].value = -1.0;
    else
      model->coeffs[i].value = 0.0;
  }

  double solution[model->num_constraints];
  SimplexLoop(model, solution);

  // Check if any artificial is still basic with positive value
  bool feasible = true;
  for (size_t i = 0; i < model->num_constraints; i++) {
    int basic = model->basics_vector[i];
    if (model->coeffs[basic].type == ARTIFICIAL && solution[i] > 1e-6) {
      feasible = false;
      break;
    }
  }

  // Restore real objective coefficients
  for (int i = 0; i < total_cols; i++)
    model->coeffs[i].value = saved[i];
  free(saved);


  for (int i = 0; i < model->artificials_count; i++) {
    int art_col = model->artificials_vector[i];
    model->coeffs[art_col].value = 0.0;
    for (int j = 0; j < model->num_constraints; j++)
      model->constraints[j].lhs_vector[art_col] = 0.0;
  }

  if (!feasible) {
    printf("Model is infeasible.\n");
    return 0;
  }
  printf("Model is feasible as of phase 1 \n");
  return 1;
}

int SimplexLoop(Model *model, double *solution_out) {




  size_t n = model->num_constraints;
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1;
  double **B_inv = Get_BasisInverse(model, model->solver_iterations);

  double original_RHS[n];
  for (size_t i = 0; i < n; i++)
    original_RHS[i] = model->constraints[i].rhs;
  // UpdateRhs(model, original_RHS, B_inv);

  while (1) {
    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      for (size_t i = 0; i < n; i++) free(B_inv[i]);
      free(B_inv);
      return -2;
    }

    printf("Beginning solver iteration %i\n", model->solver_iterations);

    double *Simplex_multiplier = Get_SimplexMultiplier(model, B_inv);

    int entering_var_idx = 0, entering_var = 0;
    double best_reduced_cost = -DBL_MAX;
    int feasibility_check = 0;

    for (size_t i = 0; i < model->non_basics_count; i++) {
      int non_basic_idx = model->non_basics[i];
      double rc = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
      if (rc > best_reduced_cost) {
        best_reduced_cost = rc;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (rc <= 0)
        feasibility_check++;
    }

    if (feasibility_check == model->non_basics_count) {
      if (solution_out)
        for (size_t i = 0; i < n; i++)
          solution_out[i] = original_RHS[i];
      for (size_t i = 0; i < n; i++) free(B_inv[i]);
      free(B_inv);
      free(Simplex_multiplier);
      return 1;
    }

    double *Pivot = Get_pivot_column(B_inv, model, entering_var);

    int exiting_var_idx = -1, exiting_var = -1;
    double best_ratio = DBL_MAX;

    for (size_t i = 0; i < n; i++) {
      if (Pivot[i] <= 1e-6) continue;
      double ratio = original_RHS[i] / Pivot[i];
      if (ratio < best_ratio && ratio >= 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }

    if (exiting_var_idx == -1) {
      free(Pivot);
      free(Simplex_multiplier);
      for (size_t i = 0; i < n; i++) free(B_inv[i]);
      free(B_inv);
      return -1;
    }

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;

    Update_BasisInverse(B_inv, Pivot, exiting_var_idx, n);

    if (model->solver_iterations % REINVERSION_FREQ == 0) {
      for (size_t i = 0; i < n; i++) free(B_inv[i]);
      free(B_inv);
      B_inv = Get_BasisInverse(model, model->solver_iterations);

      for (size_t i = 0; i < n; i++)
        original_RHS[i] = model->constraints[i].rhs;
      UpdateRhs(model, original_RHS, B_inv);
    } else {

      for (size_t i = 0; i < n; i++) {
        if (i == (size_t)exiting_var_idx)
          original_RHS[i] = best_ratio;
        else
          original_RHS[i] -= Pivot[i] * best_ratio;
      }
    }

    model->solver_iterations++;
    free(Pivot);
    free(Simplex_multiplier);
  }
}
