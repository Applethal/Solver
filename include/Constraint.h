#ifdef CONSTRAINT_H
#define CONSTRAINT_H

typedef struct Constraint {
  double *lhs_vector;
  double rhs;
  char *constraints_symbols;
  int *artificials_vector;
  int *basics_vector;
  int *non_basics;

  int basics_count;
  int non_basics_count;
  int artificials_count;
  
} Constraint ;




#endif // CONSTRAINT_H

