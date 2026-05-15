#include "core.h"
#include <float.h>




Model *ReadCsv(FILE *csvfile) {
  Model *model = (Model *)malloc(sizeof(Model));
  char line[MAX_LINE_LENGTH];

  // Line 1: Read objective type, num_vars, num_constraints
  if (!fgets(line, sizeof(line), csvfile)) {
    fprintf(stderr, "Error: Could not read header line\n");
    free(model);
    return NULL;
  }

  size_t len = strlen(line);
  if (len == sizeof(line) - 1 && line[len - 1] != '\n') {
    fprintf(stderr, "Error: Line too long (max %d chars)\n", MAX_LINE_LENGTH);
    free(model);
    return NULL;
  }

  line[strcspn(line, "\r\n")] = 0;

  char *token = strtok(line, ",");
  if (!token) {
    fprintf(stderr, "Error: Invalid header format\n");
    free(model);
    return NULL;
  }

  if (strcmp(token, "MINIMIZE") == 0) {
    model->objective = -1;
  } else if (strcmp(token, "MAXIMIZE") == 0) {
    model->objective = 1;
  } else {
    fprintf(stderr,
            "Error: Invalid objective type (expected MINIMIZE or MAXIMIZE)\n");
    free(model);
    return NULL;
  }
  token = strtok(NULL, ",");
  if (!token) {
    fprintf(stderr, "Error: Missing variables count in header\n");
    free(model);
    return NULL;
  }
  model->num_vars = atoi(token);

  token = strtok(NULL, ",");
  if (!token) {
    fprintf(stderr, "Error: Missing constraints count in header\n");
    free(model);
    return NULL;
  }
  model->num_constraints = atoi(token);

  model->bounded_vars = 0;
token = strtok(NULL, ",");
if (token) {
    model->bounded_vars = atoi(token);
}

  if (model->num_vars > MAX_VARS) {
    fprintf(stderr, "Error: Too many variables (%zu). Maximum allowed: %d\n",
            model->num_vars, MAX_VARS);
    free(model);
    return NULL;
  }

  if (model->num_constraints > MAX_CONSTRAINTS) {
    fprintf(stderr, "Error: Too many constraints (%zu). Maximum allowed: %d\n",
            model->num_constraints, MAX_CONSTRAINTS);
    free(model);
    return NULL;
  }

  // Line 2: Read objective coefficients
  if (!fgets(line, sizeof(line), csvfile)) {
    fprintf(stderr, "Error: Could not read objective coefficients\n");
    free(model);
    return NULL;
  }

  // Allocate and parse coefficients
  model->coeffs = (Variable *)malloc(model->num_vars * sizeof(Variable));
  model->integer_vars_idx = (int *)calloc(model->num_vars, sizeof(int));
  token = strtok(line, ",");
  model->integer_vars_count = 0;
  int idx = 0;
  model->bigM = 1;
  while (token != NULL && idx < model->num_vars) {
    char *i_marker = strchr(token, 'I');
    char *b_marker = strchr(token, 'B');
    if (i_marker != NULL) {
      *i_marker = '\0';
      model->coeffs[idx].value = atof(token) * model->objective;
      model->coeffs[idx].type = INTEGER;
      model->coeffs[idx].lb = 0; 
      model->coeffs[idx].ub = DBL_MAX; 
      model->integer_vars_idx[model->integer_vars_count] = idx;
      model->integer_vars_count += 1;
    } else if (b_marker != NULL) {
      *b_marker = '\0';
      model->coeffs[idx].value = atof(token) * model->objective;
      model->coeffs[idx].type = BINARY;
      model->coeffs[idx].ub = 1;
      model->coeffs[idx].lb = 0;
      model->integer_vars_idx[model->integer_vars_count] = idx;
      model->integer_vars_count += 1;
      
    } else {
      model->coeffs[idx].value = atof(token) * model->objective;
      model->coeffs[idx].type = STANDARD;
      model->coeffs[idx].lb = 0;
      model->coeffs[idx].ub = DBL_MAX;
    }

    if (model->bigM < model->coeffs[idx].value) {
      model->bigM = model->coeffs[idx].value;
    }
    idx++;

    token = strtok(NULL, ",");
  }
  model->integer_vars_idx =
      realloc(model->integer_vars_idx, model->integer_vars_count * sizeof(int));

  if (idx != model->num_vars) {
    fprintf(stderr, "Error: Expected %zu variable coefficients but got %d\n",
            model->num_vars, idx);
    free(model->coeffs);
    free(model);
    return NULL;
  }

  // Allocate arrays for raw constraint data
  model->lhs_matrix =
      (double **)malloc(model->num_constraints * sizeof(double *));
  for (int i = 0; i < model->num_constraints; i++) {
    model->lhs_matrix[i] = (double *)malloc(model->num_vars * sizeof(double));
  }

  model->rhs_vector = (double *)malloc(model->num_constraints * sizeof(double));
  model->constraints_symbols =
      (char *)malloc(model->num_constraints * sizeof(char));
  int num_le = 0, num_ge = 0, num_eq = 0;

  // Read all constraints in one pass
  for (int i = 0; i < model->num_constraints; i++) {
    if (!fgets(line, sizeof(line), csvfile)) {
      fprintf(
          stderr,
          "Error: Unexpected end of file at constraint %d. This might indicate "
          "a missing constraint or wrong constraints count in the header\n",
          i);
      for (int j = 0; j <= i; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model);
      return NULL;
    }

    token = strtok(line, ",");
    int col_idx = 0;

    while (token != NULL && col_idx < model->num_vars) {
      model->lhs_matrix[i][col_idx] = atof(token);
      col_idx++;
      token = strtok(NULL, ",");
    }

    if (col_idx != model->num_vars) {
      fprintf(
          stderr,
          "Error: Constraint %d has %d coefficients, expected %zu variables\n",
          i, col_idx, model->num_vars);
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model);
      return NULL;
    }

    if (token == NULL) {
      fprintf(stderr, "Error: Missing operator in constraint %d\n", i);
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model);
      return NULL;
    }

    char op[4] = "";
    int op_idx = 0;
    for (char *p = token; *p && op_idx < 3; p++) {
      if (*p != ' ' && *p != '\t') {
        op[op_idx++] = *p;
      }
    }
    op[op_idx] = '\0';

    if (strcmp(op, "<=") == 0) {
      model->constraints_symbols[i] = 'L';
      num_le++;
    } else if (strcmp(op, ">=") == 0) {
      model->constraints_symbols[i] = 'G';
      num_ge++;
    } else if (strcmp(op, "=") == 0) {
      model->constraints_symbols[i] = 'E';
      num_eq++;
    } else {
      fprintf(stderr, "Error: Invalid operator '%s' in constraint %d\n", op, i);
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model);
      return NULL;
    }

    token = strtok(NULL, ",");
    if (token == NULL) {
      fprintf(stderr, "Error: Missing RHS in constraint %d\n", i);
      for (int j = 0; j < model->num_constraints; j++) {
        free(model->lhs_matrix[j]);
      }
      free(model->lhs_matrix);
      free(model->rhs_vector);
      free(model->constraints_symbols);
      free(model->coeffs);
      free(model);
      return NULL;
    }
    model->rhs_vector[i] = atof(token);
  }

  // Read bounds if any were declared in the header
  if (model->bounded_vars > 0) {
    if (!fgets(line, sizeof(line), csvfile)) {
      fprintf(stderr, "Error: Expected BOUNDS line but reached end of file\n");
      // cleanup omitted for brevity — mirror your pattern above
      return NULL;
    }
    line[strcspn(line, "\r\n")] = 0;
    if (strncmp(line, "BOUNDS", 6) != 0) {
      fprintf(stderr, "Error: Expected 'BOUNDS' sentinel but got '%s'\n", line);
      return NULL;
    }

    for (int i = 0; i < model->bounded_vars; i++) {
      if (!fgets(line, sizeof(line), csvfile)) {
        fprintf(stderr, "Error: Unexpected end of file at bound %d\n", i);
        return NULL;
      }
      token = strtok(line, ",");
      if (!token) { fprintf(stderr, "Error: Missing var index in bound %d\n", i); return NULL; }
      int var_idx = atoi(token) - 1;  // 1-based in file, 0-based indexing in the program

      if (var_idx < 0 || var_idx >= model->num_vars) {
        fprintf(stderr, "Error: Bound %d references out-of-range variable %d\n", i, var_idx + 1);
        return NULL;
      }

      token = strtok(NULL, ",");
      if (!token) { fprintf(stderr, "Error: Missing lower bound for var %d\n", var_idx + 1); return NULL; }
      model->coeffs[var_idx].lb = atof(token);

      token = strtok(NULL, ",");
      if (!token) { fprintf(stderr, "Error: Missing upper bound for var %d\n", var_idx + 1); return NULL; }
      char *trimmed = token + strspn(token, " \t");
      model->coeffs[var_idx].ub = (strncmp(trimmed, "inf", 3) == 0) ? DBL_MAX : atof(token);
    }
  }

  model->artificials_vector = NULL;
  model->basics_vector = NULL;
  model->non_basics = NULL;
  model->non_basics_count = 0;
  model->solver_iterations = 1;
  model->objective_function = 0.0;
  model->slacks_surplus_count = num_le + num_ge;
  model->artificials_count = num_eq + num_ge;

  // if (model->num_constraints > model->num_vars)
  // {
  //   printf("Warning: More constraints than variables! The revised simplex
  //   algorithm works best when there are more variables than constraints!\n");
  // }
  //
return model;
}
