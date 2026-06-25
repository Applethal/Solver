#include "core.h"
#include "Variable.h"
#include "Constraint.h"
#include <stdio.h>




void TestImplementation(){

  double objective_functions[10] = {11.828861, 55.750000, 26.00000, 20.750000}; // TODO: Get all Objective functions and record them here 
  bool Debug = 0;
  double results[10] = {0};

  char *paths[] = {"../Instances/gurobi_diet.csv",
                  "../Instances/bounded/test7.csv",
                  "../Instances/bounded/test8.csv", 
                  "../Instances/bounded/small.csv"};
  

  // replace 10 with actual size of paths array
  for (int i = 0; i < 10; i++) {

    FILE *file = fopen(paths[i], "r");
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

    results[i] = model->objective_function;
    FreeModel(model);

  }
  
  for (int i = 0; i < 10; i++) {
    if (fabs(objective_functions[i] - results[i]) < 1e-4) {
      
        printf("Test passed!");

    }else {
      printf("Test case %i comes with an incorrect objective function, result = %f, expected = %f", i, results[i], objective_functions[i]);
    }

  
  }



}


