#ifndef CORE_H
#define CORE_H

#include <stddef.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "Variable.h"
#include "Constraint.h"



// Bit flags just to help decide which solving mode to go for 
#define SOLVER_INTEGER  (1 << 0)  
#define SOLVER_BOUNDED  (1 << 1)  
#define SOLVER_DEBUG    (1 << 2) // Debug support available only for BigM Solver, I have no plans yet for implementing it for the two phase method 




// I was planning to make an API out of this code logic for an online interface
// hence why I defined some limits on what the user can give to the code, feel
// free to adjust these
#define MAX_LINE_LENGTH 1024
#define MAX_VARS 2000
#define MAX_CONSTRAINTS 2000
#define MAX_MEM_BYTES 1024
#define REINVERSION_FREQ 50 // Initially used 25, 50 should be okay too



typedef struct Model {
  size_t num_constraints;      // Number of constraints
  size_t num_vars;             // Number of variables
  size_t bounded_vars;         // Number of bounded variables
  double objective_function;   
  double bigM;                 // max (coeffs) * 2
  Variable *coeffs;            // Variable objective coefficients
  Constraint *constraints;     // constraints vector 
  unsigned int *basics_vector;
  unsigned int *non_basics;
  unsigned int *artificials_vector;
  unsigned int *integer_vars_idx;
  unsigned int slacks_surplus_count;    // Slacks and surplus are mechanically the same 
  unsigned int solver_iterations;      
  unsigned int integer_vars_count;       
  unsigned int non_basics_count;
  unsigned int artificials_count;
  signed char objective;              // MINIMIZE -1, MAXIMIZE 1
} Model;

// Function declarations
void MainLoop(int argc, char* argv[]);
void PrintHelp();
Model *ReadCsv(FILE *csvfile);
void TransformModel(Model *model);
void PrintConstraints_matrix(Model *model); // Prints matrix in canonical form 
void PrintConstraints_matrix_original(Model *model); // Prints matrix prior to transformation
double **Get_BasisInverse(Model *model, int iteration);
void InvertMatrix(double **matrix, size_t n); // I no longer use this. I use rank 1 matrix inversion
void RevisedSimplex(Model *model); 
void RevisedSimplex_Debug(Model *model); 
double Get_ReducedPrice(Model *model, double **B_inv, int var_col, double *multiplier_vector);
void Get_SimplexMultiplier(Model *model, double **B_inv, double *Simplex_multiplier);
void Get_pivot_column(double **B_inv, Model *model, int best_cost_idx, double *Pivot);
void UpdateRhs(Model *model, double *rhs_vector_copy, double **B);
void Get_ObjectiveFunction(Model *model, double *rhs_vector);
void FreeModel(Model *model);
void ValidateModelPointers(Model *model);
size_t ModelMemSize(Model *model);



void TestImplementation(); // In case I want to make a series of tests, TODO: Get every proven optimal solution





// Bounded variable Simplex 

void SolveBounded(Model* model); 
void BoundedSimplex(Model* model);
void TransformBoundedModel(Model* model); // Bounded variable simplex will have a different procedure
void Get_ObjectiveFunctionBounded(Model *model, double *rhs_vector);
void Update_BasisInverse(double **B_inv, double *Pivot, int pivot_row, int n); 
void Bounded_UpdateRhs(Model *model, double *original_RHS, double *current_RHS, double **B); 


// Two-Phase code
int  SimplexLoop(Model *model, double *solution_out); 
int  RunPhase1(Model *model);   // returns 1 feasible, 0 infeasible
void RunPhase2(Model *model);   // optimizes from Phase 1 basis
void Solve(Model *model);       // Solver routine using the Two-Phase scheme
void Solve_BigM(Model *model); // Solver routine using BigM, keeping it here for the sake of legacy
void Solve_BigM_Debug(Model *model);





#endif // CORE_H
