#include "core.h"
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void PrintHelp() {
  printf("NAME\n");
  printf("     RSA - Linear programming solver using the Revised Simplex "
         "Algorithm\n\n");
  printf("SYNOPSIS\n");
  printf("     ./RSA csv file path [-Debug]\n\n");
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
  printf("     Visit the Github repo https://github.com/Applethal/RSA to see "
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
  printf("             steps during the solving process. \n");
  printf("EXIT STATUS\n");
  printf("     0       Optimal solution obtained.\n");
  printf("     1       File opening error (file does not exist or cannot be "
         "accessed).\n");
  printf("     Other   Model is infeasible or unbounded.\n\n");
  printf("OUTPUT\n");
  printf("     The program displays:\n");
  printf("     - Start message with usage hint\n");
  printf("     - Debug information (if -Debug flag is used)\n");
  printf("     - RSA iterations count\n");
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
         "https://github.com/Applethal/RSA\n\n");
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
  printf("Starting solver, to enable iterative debugging add the flag '-Debug' "
         "as an argument after the file path.\n\n");
  fclose(file);

  if (model->integer_vars_count > 0) {
    printf("Integer programming detected\n");
    IntegerSolvingLoop(model);
  } else {
    if (Debug) {
      Solve_BigM_Debug(model);
    } else {
      //Solve(model);
      Solve_BigM(model);
      //SolveBounded(model);
    }
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

void TransformModel(Model *model) {
  // TODO: Make a seperate function for two phase to make it more efficient by directly attributing bigM to the artificials for the BigM method and -1 for the two phase method
// TODO2: Add bounds to slacks
  if (!model) {
    fprintf(stderr, "Error: NULL model pointer\n");
    exit(1);
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
      model->coeffs[artif_col].value = -(model->bigM * 2 * model->num_vars);
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
      model->coeffs[artif_col].value = -(model->bigM * 2 * model->num_vars);
      model->coeffs[artif_col].constraint_idx = 0;
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
void PrintConstraints_matrix_original(Model *model) {

  printf("Constraints in original form:\n");
  for (int i = 0; i < model->num_constraints; i++) {
    printf("Constraint %d: ", i + 1);
    for (int j = 0; j < model->num_vars; j++) {
      printf("%5.1f ", model->lhs_matrix[i][j]);
    }
    printf("| RHS: %5.1f", model->rhs_vector[i]);
    printf("\n");
  }
  printf("\n");

  printf("Constraint symbols in the original problem:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    printf(" %c ", model->constraints_symbols[i]);
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
      printf("%5.1f ", model->lhs_matrix[i][j]);
    }
    printf("| RHS: %5.1f", model->rhs_vector[i]);
    printf("\n");
  }
  printf("\n");

  printf("Constraint symbols in the original problem:\n");
  for (size_t i = 0; i < model->num_constraints; i++) {
    printf(" %c ", model->constraints_symbols[i]);
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
        B[i][j] = model->lhs_matrix[i][basis_col];
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
      exit(1);
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
    dot_product += multiplier_vector[i] * model->lhs_matrix[i][var_col];
  }

  double reduced_cost = model->coeffs[var_col].value - dot_product;
  return reduced_cost;
}

double *Get_SimplexMultiplier(Model *model, double **B_inv) {
  size_t n = model->num_constraints;
  double *multiplier_vector = (double *)malloc(sizeof(double) * n);

  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      int basic_col_idx = model->basics_vector[j];
      sum += model->coeffs[basic_col_idx].value * B_inv[j][i];
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
      sum += B_inv[i][j] * model->lhs_matrix[j][best_cost_idx];
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

    model->objective_function += (model->coeffs[basic_idx].value * model->objective) * rhs_vector[i];

    // Print only standard variables with non-zero values
    if (model->coeffs[basic_idx].type == STANDARD && rhs_vector[i] > 1e-6) {
      // Convert back to original coefficient if it's a minimization problem
      // originally
      double original_coeff = model->coeffs[basic_idx].value;
      

      // TODO: Instead of having an IF statement, maybe make it so that MAX problems are set to 1 and MIN problems are set to -1 in model->objective 
      // Sadly, this will prevent any opportunity for me to use boolean memory size, I can somehow use a makeshift boolean by using a single char 
      // value of 1 and -1 since these are the only expected values of model->objective, so I should not be afraid in this area.
      // if (model->objective == 0) {
      //   original_coeff *= -1;
      // }
      //
      printf("Value of variable x%i is %f with coefficient %f\n", basic_idx,
             rhs_vector[i], original_coeff * model->objective);
    }
  }

  // // Convert objective function back for minimization
  // if (model->objective == 0) {
  //   model->objective_function *= -1;
  // }

  printf("Optimal solution found! Objective value: %f\n",
         model->objective_function);
}


void FreeModel(Model *model) {

  for (int i = 0; i < model->num_constraints; i++)
    free(model->lhs_matrix[i]);
  free(model->lhs_matrix);
  free(model->coeffs);
  free(model->rhs_vector);
  free(model->basics_vector);
  free(model->constraints_symbols);
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

  if (!model->lhs_matrix) {
    fprintf(stderr, "Fatal Error: model->lhs_matrix is NULL.\n");
    exit(1);
  }

  if (!model->coeffs) {
    fprintf(stderr, "Fatal Error: model->coeffs is NULL.\n");
    exit(1);
  }

  if (!model->rhs_vector) {
    fprintf(stderr, "Fatal Error: model->rhs_vector is NULL.\n");
    exit(1);
  }

  if (!model->basics_vector) {
    fprintf(stderr, "Fatal Error: model->basics_vector is NULL.\n");
    exit(1);
  }

  if (!model->constraints_symbols) {
    fprintf(stderr, "Fatal Error: model->constraints_symbols is NULL.\n");
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
    if (!model->lhs_matrix[i]) {
      fprintf(stderr, "Fatal Error: model->lhs_matrix[%d] is NULL.\n", i);
      exit(1);
    }
  }
}




  size_t ModelMemSize(Model * model) {

    size_t size = 0;
    size += sizeof(bool); // objective mode
    size += sizeof(int) * model->num_constraints * 2;
    size += sizeof(int) * model->num_vars;
    // variable struct
    size += sizeof(Variable) * model->num_vars; // struct
    size += sizeof(int) * model->num_vars;      // type
    size += sizeof(int) * model->num_vars;      // constraint index
    size += sizeof(double) * model->num_vars;   // value
    int total_cols = model->num_vars + model->slacks_surplus_count +
                     model->artificials_count;
    for (int i = 0; i < model->num_constraints; i++)
      size += sizeof(double) * total_cols;
    size += sizeof(Variable) * total_cols;           // Variables
    size += sizeof(double) * model->num_constraints; // RHS vector
    size += sizeof(int) * model->num_constraints;    // Basics vector
    size += sizeof(double);                          // Objective function
    size += sizeof(char) * model->num_constraints;   // Constraints symbols
    size += sizeof(int) * 2; // Inequalities and equalities count
    size += sizeof(int) * model->artificials_count;  // Artificials vector
    size += sizeof(int);                             // Solver iterations
    size += sizeof(int) * model->non_basics_count;   // Non-basics
    size += sizeof(int) * model->integer_vars_count; // Integer variable indices
    size += sizeof(double); // BigM
    size += sizeof(double); // Objective Constant



    // printf("Matrix inversion memory usage per iteration: %zu bytes (%.2f
    // KB)\n",
    //        inversion_memory, inversion_memory / 1024.0);
    // printf("Simplex multiplier vector size per iteration: %zu bytes (%.2f
    // KB)\n", constraints_vector, constraints_vector / 1024.0); printf("Memory
    // usage when updating the RHS vector per iteration: %zu bytes (%.2f KB)\n",
    // constraints_vector, constraints_vector / 1024.0); printf("Total dynamic
    // memory usage in each iteration: %zu bytes (%.2f KB)\n",
    // constraints_vector
    // * 2 + inversion_memory + size, (constraints_vector * 2 + inversion_memory
    // + size) / 1024.0); printf("Total iterations: %d\n",
    // model->solver_iterations);

    return size;
  }



  Model *deep_copy_model(const Model *model) {
    if (!model)
      return NULL;

    Model *model_copy = calloc(1, sizeof(Model));
    if (!model_copy)
      return NULL;

    size_t total_vars = model->num_vars + model->slacks_surplus_count +
                        model->artificials_count;

    model_copy->objective = model->objective;
    model_copy->num_constraints = model->num_constraints;
    model_copy->num_vars = model->num_vars;
    model_copy->objective_function = model->objective_function;
    model_copy->slacks_surplus_count = model->slacks_surplus_count;
    model_copy->artificials_count = model->artificials_count;
    model_copy->solver_iterations = model->solver_iterations;
    model_copy->non_basics_count = model->non_basics_count;
    model_copy->integer_vars_count = model->integer_vars_count;

    if (model->integer_vars_idx) {
      model_copy->integer_vars_idx =
          calloc(model->integer_vars_count, sizeof(int));
      if (!model_copy->integer_vars_idx)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->integer_vars_idx, model->integer_vars_idx,
             model_copy->integer_vars_count * sizeof(int));
    }

    if (model->lhs_matrix) {
      model_copy->lhs_matrix = calloc(model->num_constraints, sizeof(double *));
      if (!model_copy->lhs_matrix)
        return FreeModel(model_copy), NULL;

      for (size_t i = 0; i < model->num_constraints; i++) {
        model_copy->lhs_matrix[i] = calloc(total_vars, sizeof(double));
        if (!model_copy->lhs_matrix[i])
          return FreeModel(model_copy), NULL;
        memcpy(model_copy->lhs_matrix[i], model->lhs_matrix[i],
               total_vars * sizeof(double));
      }
    }

    if (model->coeffs) {
      model_copy->coeffs = calloc(total_vars, sizeof(Variable));
      if (!model_copy->coeffs)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->coeffs, model->coeffs, total_vars * sizeof(Variable));
    }

    if (model->rhs_vector) {
      model_copy->rhs_vector = calloc(model->num_constraints, sizeof(double));
      if (!model_copy->rhs_vector)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->rhs_vector, model->rhs_vector,
             model->num_constraints * sizeof(double));
    }

    if (model->basics_vector) {
      model_copy->basics_vector = calloc(model->num_constraints, sizeof(int));
      if (!model_copy->basics_vector)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->basics_vector, model->basics_vector,
             model->num_constraints * sizeof(int));
    }

    if (model->constraints_symbols) {
      model_copy->constraints_symbols =
          calloc(model->num_constraints, sizeof(char));
      if (!model_copy->constraints_symbols)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->constraints_symbols, model->constraints_symbols,
             model->num_constraints * sizeof(char));
    }

    if (model->artificials_vector && model->artificials_count > 0) {
      model_copy->artificials_vector =
          calloc(model->artificials_count, sizeof(int));
      if (!model_copy->artificials_vector)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->artificials_vector, model->artificials_vector,
             model->artificials_count * sizeof(int));
    }

    if (model->non_basics && model->non_basics_count > 0) {
      model_copy->non_basics = calloc(model->non_basics_count, sizeof(int));
      if (!model_copy->non_basics)
        return FreeModel(model_copy), NULL;
      memcpy(model_copy->non_basics, model->non_basics,
             model->non_basics_count * sizeof(int));
    }

    return model_copy;
  }


