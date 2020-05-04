SET(OPTIMIZATION_HEADERS
	case.h
	case_handler.h
	case_transfer_object.h
	#
	optimizer.h
	normalizer.h
	hybrid_optimizer.h
	optimization_exceptions.h
	#
	objective/NPV.h
	objective/objective.h
	objective/weightedsum.h
	#
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
	#
	optimizers/APPS.h
	optimizers/ExhaustiveSearch2DVert.h
	optimizers/GSS.h
	optimizers/GeneticAlgorithm.h
	optimizers/PSO.h
	optimizers/CMA_ES.h
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
	optimizers/trust_region/TrustRegionOptimization.h
	optimizers/trust_region/TrustRegionModel.h
	optimizers/trust_region/TrustRegionMath.h
	optimizers/ensemble_exp_value.h
	#
	../ThirdParty/snopt/handlers/SNOPTHandler.h
)

SET(OPTIMIZATION_SOURCES
	case.cpp
	case_handler.cpp
	case_transfer_object.cpp
	#
	optimizer.cpp
	normalizer.cpp
	hybrid_optimizer.cpp
	#
	objective/NPV.cpp
	objective/objective.cpp
	objective/weightedsum.cpp
	#
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
	#
	optimizers/APPS.cpp
	optimizers/ExhaustiveSearch2DVert.cpp
	optimizers/GSS.cpp
	optimizers/GeneticAlgorithm.cpp
	optimizers/PSO.cpp
	optimizers/CMA_ES.cpp
	optimizers/RGARDD.cpp
	optimizers/VFSA.cpp
  optimizers/SPSA.cpp
	optimizers/bayesian_optimization/AcquisitionFunction.cpp
	optimizers/bayesian_optimization/EGO.cpp
	optimizers/bayesian_optimization/af_optimizers/AFCompassSearch.cpp
	optimizers/bayesian_optimization/af_optimizers/AFOptimizer.cpp
	optimizers/bayesian_optimization/af_optimizers/AFPSO.cpp
	optimizers/compass_search.cpp
	optimizers/trust_region/TrustRegionOptimization.cpp
	optimizers/trust_region/TrustRegionModel.cpp
	optimizers/trust_region/TrustRegionMath.cpp
	optimizers/ensemble_exp_value.cpp
	#
	../ThirdParty/snopt/handlers/SNOPTHandler.cpp
	../ThirdParty/snopt/handlers/SNOPTLoader.c
	../ThirdParty/snopt/handlers/LibraryHandler.c
	solvers/SNOPTSolver.cpp
)

SET(OPTIMIZATION_TESTS
	tests/test_resource_cases.h
	tests/test_resource_optimizer.h
	tests/test_resource_test_functions.h
	#
	tests/constraints/test_bhp_constraint.cpp
	tests/constraints/test_constraint_handler.cpp
	tests/constraints/test_interwell_distance.cpp
	tests/constraints/test_pseudo_cont_boundary_2d.cpp
	tests/constraints/test_rate_constraint.cpp
	tests/constraints/test_reservoir_boundary.cpp
	tests/constraints/test_spline_well_length.cpp
	#
	tests/test_case.cpp
	tests/test_case_handler.cpp
	tests/test_case_transfer_object.cpp
	tests/test_normalizer.cpp
	tests/objective/test_weightedsum.cpp
	#
	tests/optimizers/test_apps.cpp
	tests/optimizers/test_compass_search.cpp
	tests/optimizers/test_ego.cpp
	tests/optimizers/test_ga.cpp
	tests/optimizers/test_pso.cpp
	tests/optimizers/test_vfsa.cpp
	tests/optimizers/test_spsa.cpp
	tests/optimizers/test_cma_es.cpp
	tests/optimizers/test_tr-dfo.cpp
	tests/optimizers/test_tr-dfo_exp-value.cpp
	tests/optimizers/test_tr-model-data.hpp
	tests/optimizers/test_tr-support.hpp
)
