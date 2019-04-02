SET(OPTIMIZATION_HEADERS
	case.h
	case_handler.h
	case_transfer_object.h
	constraints/bhp_constraint.h
	constraints/combined_spline_length_interwell_distance.h
	constraints/combined_spline_length_interwell_distance_reservoir_boundary.h
	constraints/constraint.h
	constraints/constraint_handler.h
	constraints/icv_constraint.h
	constraints/interwell_distance.h
	constraints/packer_constraint.h
	constraints/pseudo_cont_boundary_2d.h
	constraints/rate_constraint.h
	constraints/reservoir_boundary.h
	constraints/well_spline_constraint.h
	constraints/well_spline_length.h
	constraints/polar_well_length.h
	constraints/polar_azimuth.h
	constraints/polar_elevation.h
	constraints/polar_spline_boundary.h
	constraints/polar_xyz_boundary.h
	constraints/reservoir_boundary_toe.h
	constraints/reservoir_xyz_boundary.h
	hybrid_optimizer.h
	normalizer.h
	objective/NPV.h
	objective/objective.h
	objective/weightedsum.h
	optimization_exceptions.h
	optimizer.h
	optimizers/APPS.h
	optimizers/ExhaustiveSearch2DVert.h
	optimizers/GSS.h
	optimizers/GeneticAlgorithm.h
	optimizers/PSO.h
	optimizers/RGARDD.h
	optimizers/VFSA.h
    optimizers/SPSA.h
	optimizers/bayesian_optimization/AcquisitionFunction.h
	optimizers/bayesian_optimization/EGO.h
	optimizers/bayesian_optimization/af_optimizers/AFCompassSearch.h
	optimizers/bayesian_optimization/af_optimizers/AFOptimizer.h
	optimizers/bayesian_optimization/af_optimizers/AFPSO.h
	optimizers/compass_search.h
	optimizers/gss_patterns.hpp
)

SET(OPTIMIZATION_SOURCES
	case.cpp
	case_handler.cpp
	case_transfer_object.cpp
	constraints/bhp_constraint.cpp
	constraints/combined_spline_length_interwell_distance.cpp
	constraints/combined_spline_length_interwell_distance_reservoir_boundary.cpp
	constraints/constraint.cpp
	constraints/constraint_handler.cpp
	constraints/icv_constraint.cpp
	constraints/interwell_distance.cpp
	constraints/packer_constraint.cpp
	constraints/pseudo_cont_boundary_2d.cpp
	constraints/rate_constraint.cpp
	constraints/reservoir_boundary.cpp
	constraints/well_spline_constraint.cpp
	constraints/well_spline_length.cpp
	constraints/polar_well_length.cpp
	constraints/polar_azimuth.cpp
	constraints/polar_elevation.cpp
	constraints/polar_spline_boundary.cpp
	constraints/polar_xyz_boundary.cpp
	constraints/reservoir_boundary_toe.cpp
	constraints/reservoir_xyz_boundary.cpp
	hybrid_optimizer.cpp
	normalizer.cpp
	objective/NPV.cpp
	objective/objective.cpp
	objective/weightedsum.cpp
	optimizer.cpp
	optimizers/APPS.cpp
	optimizers/ExhaustiveSearch2DVert.cpp
	optimizers/GSS.cpp
	optimizers/GeneticAlgorithm.cpp
	optimizers/PSO.cpp
	optimizers/RGARDD.cpp
	optimizers/VFSA.cpp
    optimizers/SPSA.cpp
	optimizers/bayesian_optimization/AcquisitionFunction.cpp
	optimizers/bayesian_optimization/EGO.cpp
	optimizers/bayesian_optimization/af_optimizers/AFCompassSearch.cpp
	optimizers/bayesian_optimization/af_optimizers/AFOptimizer.cpp
	optimizers/bayesian_optimization/af_optimizers/AFPSO.cpp
	optimizers/compass_search.cpp
)

SET(OPTIMIZATION_TESTS
	tests/test_resource_cases.h
	tests/test_resource_optimizer.h
	tests/test_resource_test_functions.h
	tests/constraints/test_bhp_constraint.cpp
	tests/constraints/test_constraint_handler.cpp
	tests/constraints/test_interwell_distance.cpp
	tests/constraints/test_pseudo_cont_boundary_2d.cpp
	tests/constraints/test_rate_constraint.cpp
	tests/constraints/test_reservoir_boundary.cpp
	tests/constraints/test_spline_well_length.cpp
	tests/objective/test_weightedsum.cpp
	tests/optimizers/test_apps.cpp
	tests/optimizers/test_compass_search.cpp
	tests/optimizers/test_ego.cpp
	tests/optimizers/test_ga.cpp
	tests/optimizers/test_pso.cpp
	tests/optimizers/test_vfsa.cpp
	tests/optimizers/test_spsa.cpp
	tests/test_case.cpp
	tests/test_case_handler.cpp
	tests/test_case_transfer_object.cpp
	tests/test_normalizer.cpp
)
