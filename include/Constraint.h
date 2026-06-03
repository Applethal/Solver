#ifndef CONSTRAINT_H
#define CONSTRAINT_H

typedef struct Constraint {
  double *lhs_vector;
  double rhs;
  char constraints_symbol;  // 'L', 'G', 'E'
} Constraint;




#endif // CONSTRAINT_H

