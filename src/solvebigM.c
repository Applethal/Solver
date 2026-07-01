#include "Variable.h"
#include "core.h"
#include <stdio.h>
#include "Constraint.h"

void Solve_BigM(Model *model) {
  TransformModel(model);
  ValidateModelPointers(model);
  RevisedSimplex(model);
}



  void Solve_BigM_Debug(Model * model) {
    TransformModel(model);
    ValidateModelPointers(model);
    RevisedSimplex_Debug(model);
  }


void RevisedSimplex(Model *model) {

  int termination = 0;
  size_t n = model->num_constraints;
  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1;
  double **B_inv = Get_BasisInverse(model, model->solver_iterations);
  double *solution = (double *)malloc(n * sizeof(double));
  double *Simplex_multiplier = (double *)malloc(n * sizeof(double));
  double *Pivot = (double *)malloc(n * sizeof(double));

  for (size_t i = 0; i < n; i++) {
    solution[i] = model->constraints[i].rhs;
  }

  while (termination != 1) {

    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      break;
    }

    printf("Beginning solver iteration %i ... \n", model->solver_iterations);

    Get_SimplexMultiplier(model, B_inv, Simplex_multiplier);

    int entering_var_idx = -1;
    int entering_var = -1;

    double best_reduced_cost = -DBL_MAX;
    int feasibility_check = 0;
    for (size_t i = 0; i < model->non_basics_count; i++) {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
      if (reduced_cost > best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0) {
        feasibility_check++;
      }
    }
    if (feasibility_check == model->non_basics_count) {
      printf("Solver loop terminated!\n");
      termination++;
      Get_ObjectiveFunction(model, solution);
      break;
    }

    Get_pivot_column(B_inv, model, entering_var, Pivot);

    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = DBL_MAX;
    for (size_t i = 0; i < n; i++) {
      if (Pivot[i] <= 1e-6) continue;
      double ratio = solution[i] / Pivot[i];
      if (ratio < best_ratio && ratio >= 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }
    if (exiting_var_idx == -1) {
      printf("LP is unbounded! Terminating!\n");
      break;
    }

    Update_BasisInverse(B_inv, Pivot, exiting_var_idx, n);

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;

    if (model->solver_iterations % REINVERSION_FREQ == 0) {
      for (size_t i = 0; i < n; i++) free(B_inv[i]);
      free(B_inv);
      B_inv = Get_BasisInverse(model, model->solver_iterations);
      for (size_t i = 0; i < n; i++) solution[i] = model->constraints[i].rhs;
      UpdateRhs(model, solution, B_inv);
    } else {
      for (size_t i = 0; i < n; i++) {
        if (i == (size_t)exiting_var_idx) {
          solution[i] = best_ratio;
        } else {
          solution[i] -= Pivot[i] * best_ratio;
        }
      }
    }

    model->solver_iterations++;
  }

  for (size_t i = 0; i < n; i++) free(B_inv[i]);
  free(B_inv);
  free(solution);
  free(Simplex_multiplier);
  free(Pivot);
}



void RevisedSimplex_Debug(Model *model) {

  printf("Debug Mode: On\n");
  printf("Model after transformation:\n");
  if (model->objective == -1) {
    printf("Objective function mode: MINIMIZE\n");
  } else {
    printf("Objective function mode: MAXIMIZE");
  }
  printf("Number of variables: %zu\n", model->num_vars);
  printf("Number of constraints: %zu\n", model->num_constraints);
  printf("Objective coefficients:");
  for (int i = 0; i < model->num_vars; i++) {
    printf(" %.1f", model->coeffs[i].value);
  }
  printf("\n\n");
  PrintConstraints_matrix(model);

  int termination = 0;
  size_t n = model->num_constraints;
  printf("\n");
  printf("===================================================================\n");

  printf("Initial Non basic variable indices:\n");
  for (size_t i = 0; i < model->non_basics_count; i++) {
    printf(" %i ", model->non_basics[i]);
  }
  printf("\n");
  printf("Initial Basic variable indices:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    printf(" %i ", model->basics_vector[i]);
  }
  printf("\n");
  printf("===================================================================\n");

  int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1;

  double *solution = (double *)malloc(n * sizeof(double));
  double *Simplex_multiplier = (double *)malloc(n * sizeof(double));
  double *Pivot = (double *)malloc(n * sizeof(double));

  if (model->solver_iterations == 1) {
    printf("Solver iteration 1, B^-1 is identity matrix\n");
  } else {
    printf("Getting B^-1 matrix:\n");
  }
  double **B_inv = Get_BasisInverse(model, model->solver_iterations);

  for (size_t i = 0; i < n; i++) {
    solution[i] = model->constraints[i].rhs;
  }
  printf("Initial RHS vector (B^-1 * b):\n");
  for (size_t i = 0; i < n; i++) {
    printf(" %f ", solution[i]);
  }
  printf("\n");

  while (termination != 1) {
    int feasibility_check = 0;
    if (model->solver_iterations == MAX_ITERATIONS) {
      printf("Max iterations reached. Terminating!\n");
      break;
    }

    printf("Beginning solver iteration %i ... \n", model->solver_iterations);
    printf("===================================================================\n");
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        printf(" %f ", B_inv[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    printf("===================================================================\n");

    printf("Getting the simplex multiplier P vector:\n");
    Get_SimplexMultiplier(model, B_inv, Simplex_multiplier);
    printf("P vector elements:\n");
    for (size_t i = 0; i < n; i++) {
      printf(" %f ", Simplex_multiplier[i]);
    }
    printf("\n");
    printf("===================================================================\n");

    int entering_var_idx = 0;
    int entering_var = 0;

    printf("Current RHS vector:\n");
    for (size_t i = 0; i < n; i++) {
      printf(" %f ", solution[i]);
    }
    printf("\n");

    double best_reduced_cost = -DBL_MAX;
    printf("MAXIMIZATION problem (after conversion), choosing the highest reduced cost value\n");
    printf("===================================================================\n");
    printf("Getting the entering variable:\n");
    printf("===================================================================\n");

    for (size_t i = 0; i < model->non_basics_count; i++) {
      int non_basic_idx = model->non_basics[i];
      double reduced_cost = Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
      printf("The reduced cost of non basic variable %i is %f \n", non_basic_idx, reduced_cost);

      if (reduced_cost > best_reduced_cost) {
        best_reduced_cost = reduced_cost;
        entering_var_idx = i;
        entering_var = non_basic_idx;
      }
      if (reduced_cost <= 0) {
        feasibility_check++;
      }
    }

    printf("Entering variable index: %i, cost: %f \n", entering_var, best_reduced_cost);
    printf("===================================================================\n");

    if (feasibility_check == model->non_basics_count) {
      printf("Solver loop terminated!\n");
      termination++;
      Get_ObjectiveFunction(model, solution);
      break;
    }

    Get_pivot_column(B_inv, model, entering_var, Pivot);
    printf("Getting the exiting variable:\n");
    printf("===================================================================\n");

    int exiting_var_idx = -1;
    int exiting_var = -1;
    double best_ratio = DBL_MAX;

    for (size_t i = 0; i < n; i++) {
      if (Pivot[i] <= 1e-6) continue;
      double ratio = solution[i] / Pivot[i];
      printf("Ratio of variable %i is %f which has a pivot value of %f \n",
             model->basics_vector[i], ratio, Pivot[i]);
      if (ratio < best_ratio && ratio >= 0) {
        best_ratio = ratio;
        exiting_var_idx = i;
        exiting_var = model->basics_vector[i];
      }
    }

    if (exiting_var_idx == -1) {
      printf("No exiting variable found, LP is unbounded! Terminating!\n");
      break;
    }

    printf("The exiting variable is %i with ratio of %f \n", exiting_var, best_ratio);
    printf("===================================================================\n");

    Update_BasisInverse(B_inv, Pivot, exiting_var_idx, n);

    model->non_basics[entering_var_idx] = exiting_var;
    model->basics_vector[exiting_var_idx] = entering_var;

    printf("Entering variable:%i at position %i \n", entering_var, entering_var_idx);
    printf("Exiting variable:%i at position: %i \n", exiting_var, exiting_var_idx);
    printf("Non basics updated:\n");
    for (size_t i = 0; i < model->non_basics_count; i++) {
      printf(" %i ", model->non_basics[i]);
    }
    printf("\n");
    printf("Basics vector:\n");
    for (size_t i = 0; i < model->num_constraints; i++) {
      printf(" %i ", model->basics_vector[i]);
    }
    printf("\n");

    if (model->solver_iterations % REINVERSION_FREQ == 0) {
      printf("\nReinversion frequency hit, rebuilding B^-1 from scratch:\n");
      for (size_t i = 0; i < n; i++) free(B_inv[i]);
      free(B_inv);
      B_inv = Get_BasisInverse(model, model->solver_iterations);
      printf("Resyncing RHS vector exactly (drift correction):\n");
      for (size_t i = 0; i < n; i++) solution[i] = model->constraints[i].rhs;
      UpdateRhs(model, solution, B_inv);
    } else {
      for (size_t i = 0; i < n; i++) {
        if (i == (size_t)exiting_var_idx) {
          solution[i] = best_ratio;
        } else {
          solution[i] -= Pivot[i] * best_ratio;
        }
      }
    }

    printf("RHS vector after update:\n");
    for (size_t i = 0; i < n; i++) {
      printf(" %f ", solution[i]);
    }
    printf("\n");

    model->solver_iterations++;

    printf("\n");
    printf("Iteration %i ends here, press any key to continue!\n", model->solver_iterations - 1);
    printf("===================================================================\n");
    getchar();
  }

  for (size_t i = 0; i < n; i++) free(B_inv[i]);
  free(B_inv);
  free(solution);
  free(Simplex_multiplier);
  free(Pivot);
}
