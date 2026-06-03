// #include "core.h"
//
//
//
//
// void IntegerSolvingLoop(Model *model) {
//
//   Model *model_copy = deep_copy_model(model);
//   TransformModel(model_copy);
//   PrintConstraints_matrix_original(model_copy);
//
//   TransformModel(model);
//   ValidateModelPointers(model);
//   int var_test = RevisedSimplex_Integer(model, NULL, NULL, NULL, NULL);
//   printf("The first non integer variable is %i \n", var_test);
//   // ModifyConstraint_IntegerUB(model_copy,
//   //                            model_copy->coeffs[var_test].constraint_idx, 10);
//
//   printf("Model after solving root node relaxation:");
//
//   PrintConstraints_matrix_original(model_copy);
//
//   TransformModel(model_copy);
//   PrintConstraints_matrix(model_copy);
//
//   exit(0);
// }
//
//
//
//
// void ModifyConstraint_IntegerUB(Model *model, int index, double value) {
//
//   model->rhs_vector[index] = value;
//
//   model->constraints_symbols[index] = 'L';
// }
// void ModifyConstraint_IntegerLB(Model *model, int index, double value) {
//
//   model->rhs_vector[index] = value;
//
//   model->constraints_symbols[index] = 'G';
// }
// void ModifyConstraint_BinaryUB(Model *model, int index) {
//
//   model->rhs_vector[index] = 1;
//
//   model->constraints_symbols[index] = 'E';
// }
// void ModifyConstraint_BinaryLB(Model *model, int index) {
//
//   model->rhs_vector[index] = 0;
//
//   model->constraints_symbols[index] = 'E';
// }
//
// int RevisedSimplex_Integer(Model *model, bool warmstart, int *parent_basis,
//                            int *parent_non_basics, double *solution_out) {
//   int termination = 0;
//   size_t n = model->num_constraints;
//   int MAX_ITERATIONS = (model->num_vars * model->num_constraints) + 1;
//
//   while (termination != 1) {
//     if (model->solver_iterations == MAX_ITERATIONS) {
//       printf("Max iterations reached. Terminating!\n");
//       FreeModel(model);
//       exit(0);
//     }
//
//     int feasibility_check = 0;
//
//     double **B_inv = Get_BasisInverse(model, model->solver_iterations);
//
//     double original_RHS[n];
//     for (size_t i = 0; i < n; i++) {
//       original_RHS[i] = model->rhs_vector[i];
//     }
//
//     double *Simplex_multiplier = Get_SimplexMultiplier(model, B_inv);
//
//     int entering_var_idx = 0;
//     int entering_var = 0;
//     UpdateRhs(model, original_RHS, B_inv);
//
//     double best_reduced_cost = -DBL_MAX;
//
//     for (size_t i = 0; i < model->non_basics_count; i++) {
//       int non_basic_idx = model->non_basics[i];
//       double reduced_cost =
//           Get_ReducedPrice(model, B_inv, non_basic_idx, Simplex_multiplier);
//
//       if (reduced_cost > best_reduced_cost) {
//         best_reduced_cost = reduced_cost;
//         entering_var_idx = i;
//         entering_var = non_basic_idx;
//       }
//       if (reduced_cost <= 0) {
//         feasibility_check++;
//       }
//     }
//
//     if (feasibility_check == model->non_basics_count) {
//       printf("Solver loop terminated!\n");
//       termination++;
//
//       for (size_t i = 0; i < n; i++) {
//         free(B_inv[i]);
//       }
//       free(B_inv);
//       free(Simplex_multiplier);
//       return CheckIntegrity(model, original_RHS);
//     }
//
//     double *Pivot = Get_pivot_column(B_inv, model, entering_var);
//
//     int exiting_var_idx = -1;
//     int exiting_var = -1;
//     double best_ratio = DBL_MAX;
//
//     for (size_t i = 0; i < n; i++) {
//       if (Pivot[i] <= 1e-6) {
//         continue;
//       }
//       double ratio = original_RHS[i] / Pivot[i];
//       if (ratio < best_ratio && ratio >= 0) {
//         best_ratio = ratio;
//         exiting_var_idx = i;
//         exiting_var = model->basics_vector[i];
//       }
//     }
//
//     if (exiting_var_idx == -1) {
//       printf("LP is unbounded! Terminating!\n");
//       exit(0);
//     }
//
//     model->non_basics[entering_var_idx] = exiting_var;
//     model->basics_vector[exiting_var_idx] = entering_var;
//     model->solver_iterations++;
//
//     free(Pivot);
//     for (size_t i = 0; i < n; i++) {
//       free(B_inv[i]);
//     }
//     free(B_inv);
//     free(Simplex_multiplier);
//   }
// }

