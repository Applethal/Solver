#ifndef VARIABLE_H
#define VARIABLE_H


typedef enum {
  STANDARD = 1,
  ARTIFICIAL = 2,
  SLACK = 3, 
  SURPLUS = 4,
  INTEGER = 5,
  BINARY = 6
} VariableType;

typedef enum { 
  LOWER, 
  UPPER, 
  BASIS
} VarStatus;




typedef struct Variable{
  double value;
  VariableType type;  
  double lb; 
  double ub;
  double originallb; 
  VarStatus status;
  double original_value; // During the bounded simplex, it might be that I will flip the signs of coefficients. I have no efficient way to track the original state for now, I hope 8 bytes won't be that much :( 

                      
} Variable;


#endif // VARIABLE_H
