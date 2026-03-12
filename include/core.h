#ifndef CORE_H
#define CORE_H

#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

typedef enum {
  STANDARD = 1,
  ARTIFICIAL = 2,
  SLACK = 3, 
  INTEGER = 4
} VariableType;

typedef struct {
  double value;
  VariableType type;  
  int constraint_idx; // -1 by default if the variable has no integrity constraints. This is going to help with the branch and bound to help track which bound constraint this variable corresponds to. For Integer variables, the value should be constraint_idx - 1 - integer variable count.
                      
} Variable;


// typedef enum {
//   INTEGER = 0,
//   INFEASIBLE = 1,
//   BRANCH = 2
// } SolverStatus;

typedef struct Subproblem{
  

  int id;
  int variable; // contains the idx of the variable to be fixed
  int value; // for debugging purposes, to track if I am enforcing the upper bound or lower bound constraint
 // SolverStatus status; // 0 = integer solution, 2 = branch, 1 = infeasible, for debugging purposes

} Subproblem;





typedef struct
{
  bool objective;           // MINIMIZE 0, MAXIMIZE 1
  size_t num_constraints;    // Number of constraints
  size_t num_vars;           // Number of variables
  double **lhs_matrix;       // Constraints Left hand side
  Variable *coeffs;          // Variable objective coefficients
  double *rhs_vector;        // Right hand side vector of the constraints
  int *basics_vector;        // Vector where I keep track of the basic variables
  double objective_function; // Self explanatory, stores the model's objective function
  char *constraints_symbols; // Tracks the constraints' symbols, for debugging purposes only
  int slacks_surplus_count;  // Counts the number of slack and surplus vars
  int artificials_count;     // Counts the number of artificial vars
  int *artificials_vector;   // Contains the indices of artificial vars
  int solver_iterations;     // Self explanatory, tracks the solver iterations count 
  int *non_basics;           // Contains indices of non-basic variables
  int non_basics_count;      // Counts the number of non-basic variables
  int integer_vars_count; // Counts the number of Integer variables. If the count is 0, then no Integer model is detected of course.
} Model;

// Function declarations
void SolverLoop(int argc, char* argv[]);
void PrintHelp();
Model *ReadCsv(FILE *csvfile);
void TransformModel(Model *model);
void Printlhs_matrix(Model *model);
double **Get_BasisInverse(Model *model, int iteration);
void InvertMatrix(double **matrix, size_t n);
void RevisedSimplex(Model *model);
void RevisedSimplex_Debug(Model *model);
double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector);
double *Get_SimplexMultiplier(Model *model, double **B_inv);
double *Get_pivot_column(double **B_inv, Model *model, int best_cost_idx);
void UpdateRhs(Model *model, double *rhs_vector_copy, double **B);
void Get_ObjectiveFunction(Model *model, double *rhs_vector);
void FreeModel(Model *model);
void ValidateModelPointers(Model *model);

#endif // CORE_H
