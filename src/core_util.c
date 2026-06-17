#include "core.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Variable.h"
#include "Constraint.h"

void PrintHelp() {
  printf("NAME\n");
  printf("     RSA - Linear programming solver using the Revised Simplex "
         "Algorithm\n\n");
  printf("SYNOPSIS\n");
  printf("     ./Solver csv file path [-Debug]\n\n");
  printf("DESCRIPTION\n");
  printf(
      "     RSA is a linear programming solver that implements the Revised\n");
  printf("     Simplex Algorithm. It reads a linear programming model from a "
         "CSV file\n");
  printf("     and computes the optimal solution for the given objective "
         "function and\n");
  printf("     constraints.\n\n");
  printf("     The program supports both maximization and minimization "
         "objectives and\n");
  printf("     automatically handles equality and inequality constraints by "
         "adding slack\n");
  printf("     and artificial variables as needed.\n\n");
  printf("     Visit the Github repo https://github.com/Applethal/Solver to see "
         "an example of how the CSV data input should look like.\n");
  printf("OPTIONS\n");
  printf("     csv file path\n");
  printf("             Path to the input CSV file containing the linear "
         "programming\n");
  printf("             model. This argument is required. The CSV file should "
         "be properly\n");
  printf("             formatted with the objective function, variables, and "
         "constraints.\n\n");
  printf("     -Debug  Enables debug mode. When this flag is provided, the "
         "program\n");
  printf("             displays detailed information about the model and shows "
         "iterative\n");
  printf("             steps during the solving process (bigM mode only). \n");
  printf("EXIT STATUS\n");
  printf("     0       Optimal solution obtained.\n");
  printf("     1       File opening error (file does not exist or cannot be "
         "accessed).\n");
  printf("     Other   Model is infeasible or unbounded.\n\n");
  printf("OUTPUT\n");
  printf("     The program displays:\n");
  printf("     - Start message with usage hint\n");
  printf("     - Debug information (if -Debug flag is used)\n");
  printf("     - Solver iterations count\n");
  printf("     - Final objective function value and the values for each "
         "variables in the final basis.\n");
  printf("     - Solving time in seconds\n");
  printf("     - Model's size in bytes\n");
  printf("AUTHOR\n");
  printf("     Written by Applethal / Saad Nordine\n\n");
  printf("     Thanks to Victor Hebert for implementing the Make file routine "
         "sequence");
  printf("REPORTING BUGS\n");
  printf("     Please report any bugs in the issues page in "
         "https://github.com/Applethal/Solver\n\n");
}

void MainLoop(int argc, char *argv[]) {
  clock_t begin = clock();
  bool Debug = false;
  if (argc < 2) {
    PrintHelp();
    exit(0);
  }
  if (argc >= 3 && strcmp(argv[2], "-Debug") == 0) {
    Debug = true;
  }
  FILE *file = fopen(argv[1], "r");
  if (file == NULL) {
    printf("Cannot open file\n");
    exit(1);
  }
  Model *model = ReadCsv(file);
  fclose(file);
  int solver_flags = 0;
  if (model->integer_vars_count > 0)
    solver_flags |= SOLVER_INTEGER;
  if (model->bounded_vars > 0)
    solver_flags |= SOLVER_BOUNDED;
  if (Debug)
    solver_flags |= SOLVER_DEBUG;



    
  switch (solver_flags & ~SOLVER_DEBUG) {
  case SOLVER_INTEGER:
    printf("Integer programming detected, still not implemented\n");
    //IntegerSolvingLoop(model);
    break;
  case SOLVER_BOUNDED:
    printf("Bounded simplex detected\n");
    if (solver_flags & SOLVER_DEBUG)
      printf("Warning: Debug mode not available for bounded simplex, running "
             "normally\n");
    SolveBounded(model);
    break;
  default:
    if (solver_flags & SOLVER_DEBUG) {
      
      Solve_BigM_Debug(model);
    } else {
      printf("Starting solver, to enable iterative debugging add the flag '-Debug' "
         "as an argument after the file path.\n\n");

      //Solve(model);
      Solve_BigM(model);
    }
    break;
  }
  printf("\n");
  printf("Objective function: %f\n", model->objective_function);
  clock_t end = clock();
  printf("Solving time: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
  size_t size = ModelMemSize(model);
  printf("Model's estimated memory usage: %zu bytes (%.2f KB)\n", size,
         size / 1024.0);
  FreeModel(model);
}

void PrintConstraints_matrix_original(Model *model) {

  printf("Constraints in original form:\n");
  for (int i = 0; i < model->num_constraints; i++) {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < model->num_vars; j++) {
      printf("%5.1f ", model->constraints[i].lhs_vector[j]);
    }
    printf("| RHS: %5.1f", model->constraints[i].rhs);
    printf("\n");
  }
  printf("\n");

  printf("Constraint symbols in the original problem:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    printf(" %c ", model->constraints[i].constraints_symbol);
  }
  printf("\n");
}
void PrintConstraints_matrix(Model *model) {
  int total_cols =
      model->num_vars + model->slacks_surplus_count + model->artificials_count;

  printf("Constraints in canonical form:\n");
  for (int i = 0; i < model->num_constraints; i++) {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < total_cols; j++) {
      printf("%5.1f ", model->constraints[i].lhs_vector[j]);
    }
    printf("| RHS: %5.1f", model->constraints[i].rhs);
    printf("\n");
  }
  printf("\n");

  printf("Constraint symbols in the original problem:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    printf(" %c ", model->constraints[i].constraints_symbol);
  }
  printf("\n");
}

double **Get_BasisInverse(Model *model, int iteration) {
  size_t n = model->num_constraints;
  double **B_inv = (double **)malloc(n * sizeof(double *));

  for (size_t i = 0; i < n; i++) {
    B_inv[i] = (double *)malloc(n * sizeof(double));
  }

  if (iteration == 1) {
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        B_inv[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }
  } else {
    double **B = (double **)malloc(n * sizeof(double *));
    for (size_t i = 0; i < n; i++) {
      B[i] = (double *)malloc(n * sizeof(double));
      for (size_t j = 0; j < n; j++) {
        int basis_col = model->basics_vector[j];
        // B[i][j] = model->constraints[i].lhs_vector[basis_col];
        B[i][j] = model->constraints[i].lhs_vector[basis_col] * model->coeffs[basis_col].flipped;

      }
    }

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        B_inv[i][j] = B[i][j];
      }
    }

    for (size_t i = 0; i < n; i++) {
      free(B[i]);
    }
    free(B);

    InvertMatrix(B_inv, n);
  }

  return B_inv;
}
//I don't need this anymore, it stays for legacy purposes
void InvertMatrix(double **matrix, size_t n) {
  const double EPS = 1e-12;

  double *block = malloc(n * (2 * n) * sizeof(double));
  double **aug = malloc(n * sizeof(double *));
  if (!block || !aug) {
    printf("ERROR: Out of memory\n");
    exit(1);
  }

  for (size_t i = 0; i < n; i++)
    aug[i] = block + i * (2 * n);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++)
      aug[i][j] = matrix[i][j];
    for (size_t j = n; j < 2 * n; j++)
      aug[i][j] = (j - n == i);
  }

  for (size_t i = 0; i < n; i++) {
    size_t max_row = i;
    double max_val = fabs(aug[i][i]);

    for (size_t k = i + 1; k < n; k++) {
      double v = fabs(aug[k][i]);
      if (v > max_val) {
        max_val = v;
        max_row = k;
      }
    }

    if (max_row != i) {
      double *tmp = aug[i];
      aug[i] = aug[max_row];
      aug[max_row] = tmp;
    }

    double pivot = aug[i][i];
    if (fabs(pivot) < EPS) {
      printf("ERROR: Singular matrix - basis is not invertible!\n");
      free(aug);
      free(block);
      exit(1); //TODO: Free memory blocks once this logic hits
    }

    double inv_pivot = 1.0 / pivot;
    for (size_t j = 0; j < 2 * n; j++)
      aug[i][j] *= inv_pivot;

    for (size_t k = 0; k < n; k++) {
      if (k == i)
        continue;
      double factor = aug[k][i];
      if (factor == 0.0)
        continue;

      for (size_t j = 0; j < 2 * n; j++)
        aug[k][j] -= factor * aug[i][j];
    }
  }

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      matrix[i][j] = aug[i][n + j];

  free(aug);
  free(block);
}
double Get_ReducedPrice(Model *model, double **B_inv, int var_col,
                        double *multiplier_vector) {
  size_t n = model->num_constraints;
  double dot_product = 0.0;

  for (int i = 0; i < n; i++) {
    dot_product += multiplier_vector[i] * model->constraints[i].lhs_vector[var_col] * model->coeffs[var_col].flipped;
  }

  // double reduced_cost = model->coeffs[var_col].value - dot_product;
  double reduced_cost = model->coeffs[var_col].value * model->coeffs[var_col].flipped - dot_product;
  return reduced_cost;
}

double *Get_SimplexMultiplier(Model *model, double **B_inv) {
  size_t n = model->num_constraints;
  double *multiplier_vector = (double *)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      int basic_col_idx = model->basics_vector[j];
      // sum += model->coeffs[basic_col_idx].value * B_inv[j][i];
      sum += model->coeffs[basic_col_idx].value * model->coeffs[basic_col_idx].flipped * B_inv[j][i];
    }
    multiplier_vector[i] = sum;
  }

  return multiplier_vector;
}

double *Get_pivot_column(double **B_inv, Model *model, int best_cost_idx) {
  size_t n = model->num_constraints;
  double *Pivot = (double *)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      // sum += B_inv[i][j] * model->constraints[j].lhs_vector[best_cost_idx];
      sum += B_inv[i][j] * model->constraints[j].lhs_vector[best_cost_idx] * model->coeffs[best_cost_idx].flipped;
    }
    Pivot[i] = sum;
  }

  return Pivot;
}

void UpdateRhs(Model *model, double *rhs_vector_copy, double **B) {
  size_t n = model->num_constraints;
  double *temp = (double *)malloc(n * sizeof(double));

  for (int i = 0; i < n; i++) {
    temp[i] = 0.0;
    for (int j = 0; j < n; j++) {
      temp[i] += B[i][j] * rhs_vector_copy[j];
    }
  }

  for (int i = 0; i < n; i++) {
    rhs_vector_copy[i] = temp[i];
  }

  free(temp);
}

int CheckIntegrity(Model *model, double *rhs_vector) {

  size_t n = model->num_constraints;

  for (int i = 0; i < n; i++) {
    int basic_idx = model->basics_vector[i];

    if (model->coeffs[basic_idx].type == ARTIFICIAL && rhs_vector[i] > 1e-6) {
      printf("Model Infeasible. Artificial basic variable has a positive RHS "
             "value. Terminating!\n");
      exit(0);
    }

    if (model->coeffs[basic_idx].type == INTEGER && rhs_vector[i] > 1e-6) {
      double original_coeff = model->coeffs[basic_idx].value * model->objective;
      // if (model->objective == 0) {
      //   original_coeff *= -1;
      // }
      if (fabs(rhs_vector[i] - floor(rhs_vector[i])) > 1e-6) {
        return basic_idx;
      }
    }

    if (model->coeffs[basic_idx].type == BINARY && rhs_vector[i] > 1e-6) {
      double original_coeff = model->coeffs[basic_idx].value * model->objective;
      // if (model->objective == 0) {
      //   original_coeff *= -1;
      // }
      if (fabs(rhs_vector[i] - floor(rhs_vector[i])) > 1e-6) {
        return basic_idx;
      }
    }
  }
  return -1;
}

void Get_ObjectiveFunction(Model *model, double *rhs_vector) {
  size_t n = model->num_constraints;

  for (int i = 0; i < n; i++) {
    int basic_idx = model->basics_vector[i];

    if (model->coeffs[basic_idx].type == ARTIFICIAL && rhs_vector[i] > 1e-6) {
      printf("Model Infeasible. Artificial basic variable has a positive RHS "
             "value. Terminating!\n");
      exit(0);
    }

    model->objective_function +=
        (model->coeffs[basic_idx].value * model->objective) * rhs_vector[i];

    // Print only standard variables with non-zero values
    if (model->coeffs[basic_idx].type == STANDARD && rhs_vector[i] > 1e-6) {
      // Convert back to original coefficient if it's a minimization problem
      // originally
      double original_coeff = model->coeffs[basic_idx].value;

      
      printf("Value of variable x%i is %f with coefficient %f\n", basic_idx,
             rhs_vector[i], original_coeff * model->objective);
    }
  }


  printf("Optimal solution found! Objective value: %f\n",
         model->objective_function);
}


void SwitchConstraint(Model *model, int index){

  for (size_t  i = 0; i < model->num_vars; i++) {
      model->constraints[index].lhs_vector[i] *= -1;

  
  }
  
  model->constraints[index].constraints_symbol = 'L';

  model->constraints[index].rhs *= -1;


}


void FreeModel(Model *model) {
  for (int i = 0; i < model->num_constraints; i++)
    free(model->constraints[i].lhs_vector);
  free(model->constraints);
  free(model->coeffs);
  free(model->basics_vector);
  free(model->artificials_vector);
  free(model->non_basics);
  free(model->integer_vars_idx);
  free(model);
  printf("-------------------------------\n");
  printf("Model free'd from the heap\n");
}

void ValidateModelPointers(Model *model) {
  if (!model) {
    fprintf(stderr, "Fatal Error: Model pointer is NULL.\n");
    exit(1);
  }
  if (!model->constraints) {
    fprintf(stderr, "Fatal Error: model->constraints is NULL.\n");
    exit(1);
  }
  if (!model->coeffs) {
    fprintf(stderr, "Fatal Error: model->coeffs is NULL.\n");
    exit(1);
  }
  if (!model->basics_vector) {
    fprintf(stderr, "Fatal Error: model->basics_vector is NULL.\n");
    exit(1);
  }
  if (!model->artificials_vector) {
    fprintf(stderr, "Fatal Error: model->artificials_vector is NULL.\n");
    exit(1);
  }
  if (!model->non_basics) {
    fprintf(stderr, "Fatal Error: model->non_basics is NULL.\n");
    exit(1);
  }
  for (int i = 0; i < model->num_constraints; i++) {
    if (!model->constraints[i].lhs_vector) {
      fprintf(stderr, "Fatal Error: model->constraints[%d].lhs_vector is NULL.\n", i);
      exit(1);
    }
  }
}

size_t ModelMemSize(Model *model) {
  size_t size = 0;
  size += sizeof(char);                            // objective mode
  size += sizeof(int) * model->num_constraints * 2;
  size += sizeof(int) * model->num_vars;
  size += sizeof(Variable) * model->num_vars;      // original coeffs struct
  size += sizeof(int) * model->num_vars;           // type
  size += sizeof(int) * model->num_vars;           // constraint index
  size += sizeof(double) * model->num_vars;        // value

  int total_cols = model->num_vars + model->slacks_surplus_count + model->artificials_count;

  // Constraint array: struct overhead + each lhs_vector
  size += sizeof(Constraint) * model->num_constraints;
  for (int i = 0; i < model->num_constraints; i++)
    size += sizeof(double) * total_cols;

  size += sizeof(Variable) * total_cols;           // expanded coeffs
  size += sizeof(int) * model->num_constraints;    // basics vector
  size += sizeof(double);                          // objective function
  size += sizeof(int) * 2;                         // inequalities and equalities count
  size += sizeof(int) * model->artificials_count;  // artificials vector
  size += sizeof(int);                             // solver iterations
  size += sizeof(int) * model->non_basics_count;   // non-basics
  size += sizeof(int) * model->integer_vars_count; // integer variable indices
  size += sizeof(double);                          // bigM
  size += sizeof(double);                          // objective constant
  return size;
}


void Update_BasisInverse(double **B_inv, double *Pivot, int pivot_row, int n) {
    // Scale the pivot row
    double pivot_elem = Pivot[pivot_row];
    for (int j = 0; j < n; j++)
        B_inv[pivot_row][j] /= pivot_elem;

    // Update all other rows: old_row -= (Pivot[i] / pivot_elem) * pivot_row
    for (int i = 0; i < n; i++) {
        if (i == pivot_row) continue;
        double factor = Pivot[i];
        for (int j = 0; j < n; j++)
            B_inv[i][j] -= factor * B_inv[pivot_row][j];
    }
}



