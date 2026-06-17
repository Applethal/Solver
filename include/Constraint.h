#ifndef CONSTRAINT_H
#define CONSTRAINT_H

typedef struct Constraint {
  double *lhs_vector;
  double rhs;
  double lhs_sum; // usable only for bounded-variable simplex. It is used to give slacks and surplus vars an upper bound
  char constraints_symbol;  // 'L', 'G', 'E'

} Constraint;




#endif // CONSTRAINT_H

