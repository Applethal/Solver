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
  double lb; 
  double ub;
  double originallb; 
  signed char flipped; // 1 if even flips, -1 if odd flips. This is my way of tracking bound flips
  VariableType type;  
  VarStatus status;
                      
} Variable;


#endif // VARIABLE_H
