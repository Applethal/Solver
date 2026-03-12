#ifndef BNB_H 
#define BNB_H 
#include <RSA.h>



typedef enum {
  INTEGER = 0,
  INFEASIBLE = 1,
  BRANCH = 2
} Status;



typedef struct Subproblem{

  int id;
  int variable; // contains the idx of the variable to be fixed
  Status status; // 1 = integer solution, 0 = branch, 2 = infeasible

} Subproblem;


void BNB(Model* model, Subproblem* subproblem);


void SubproblemCheck(Model* model,double* original_RHS,Subproblem* subproblem);




#endif // !BNB_H 
