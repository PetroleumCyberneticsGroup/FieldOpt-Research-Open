SET(OPTIMIZATION_HEADERS
		#
		# objective
		objective/NPV.h
		objective/objective.h
		objective/weightedsum.h
		#
		# case
		case.h
		case_handler.h
		case_transfer_object.h
		#
		# optimizer / optimizers
		optimizer.h
		normalizer.h
		optimization_exceptions.h
		#
		# constraints
		constraints/constraint.h
		constraints/bhp_constraint.h
		constraints/rate_constraint.h
		constraints/combined_spline_length_interwell_distance.h
		constraints/combined_spline_length_interwell_distance_reservoir_boundary.h
		constraints/reservoir_boundary.h
		constraints/interwell_distance.h
		constraints/well_spline_length.h
		constraints/well_spline_constraint.h
		#
		constraints/icv_constraint.h
		constraints/packer_constraint.h
		constraints/pseudo_cont_boundary_2d.h
		#
		constraints/polar_well_length.h
		constraints/polar_azimuth.h
		constraints/polar_elevation.h
		constraints/polar_spline_boundary.h
		constraints/polar_xyz_boundary.h
		constraints/reservoir_boundary_toe.h
		constraints/reservoir_xyz_boundary.h
		#
		# constraint-handling
		constraints/constraint_handler.h
		#
		# algorithms: pattern, stochastic
		optimizers/compass_search.h
		optimizers/APPS.h
		optimizers/gss_patterns.hpp
		optimizers/ExhaustiveSearch2DVert.h
		optimizers/GSS.h
		optimizers/GeneticAlgorithm.h
		optimizers/PSO.h
		optimizers/CMA_ES.h
		optimizers/RGARDD.h
		optimizers/VFSA.h
		optimizers/SPSA.h
		#
		# algorithms: bayesian
		optimizers/bayesian_optimization/AcquisitionFunction.h
		optimizers/bayesian_optimization/EGO.h
		optimizers/bayesian_optimization/af_optimizers/AFCompassSearch.h
		optimizers/bayesian_optimization/af_optimizers/AFOptimizer.h
		optimizers/bayesian_optimization/af_optimizers/AFPSO.h
		#
		# algorithms: df-tr
		optimizers/trust_region/TrustRegionOptimization.h
		optimizers/trust_region/TrustRegionModel.h
		optimizers/trust_region/TrustRegionMath.h
		optimizers/ensemble_exp_value.h
		#
		# algorithms: hybrid
		hybrid_optimizer.h
		#
		# solvers
		../ThirdParty/snopt/handlers/SNOPTHandler.h
		)

SET(OPTIMIZATION_SOURCES
		#
		# objective
		objective/NPV.cpp
		objective/objective.cpp
		objective/weightedsum.cpp
		#
		# case
		case.cpp
		case_handler.cpp
		case_transfer_object.cpp
		#
		# optimizer / optimizers
		optimizer.cpp
		normalizer.cpp
		#
		# constraints
		constraints/constraint.cpp
		constraints/bhp_constraint.cpp
		constraints/rate_constraint.cpp
		constraints/combined_spline_length_interwell_distance.cpp
		constraints/combined_spline_length_interwell_distance_reservoir_boundary.cpp
		constraints/reservoir_boundary.cpp
		constraints/interwell_distance.cpp
		constraints/well_spline_length.cpp
		constraints/well_spline_constraint.cpp
		#
		constraints/icv_constraint.cpp
		constraints/packer_constraint.cpp
		constraints/pseudo_cont_boundary_2d.cpp
		#
		constraints/polar_well_length.cpp
		constraints/polar_azimuth.cpp
		constraints/polar_elevation.cpp
		constraints/polar_spline_boundary.cpp
		constraints/polar_xyz_boundary.cpp
		constraints/reservoir_boundary_toe.cpp
		constraints/reservoir_xyz_boundary.cpp
		#
		# constraint-handling
		constraints/constraint_handler.cpp
		#
		# algorithms: pattern, stochastic
		optimizers/compass_search.cpp
		optimizers/APPS.cpp
		optimizers/ExhaustiveSearch2DVert.cpp
		optimizers/GSS.cpp
		optimizers/GeneticAlgorithm.cpp
		optimizers/PSO.cpp
		optimizers/CMA_ES.cpp
		optimizers/RGARDD.cpp
		optimizers/VFSA.cpp
		optimizers/SPSA.cpp
		#
		# algorithms: bayesian
		optimizers/bayesian_optimization/AcquisitionFunction.cpp
		optimizers/bayesian_optimization/EGO.cpp
		optimizers/bayesian_optimization/af_optimizers/AFCompassSearch.cpp
		optimizers/bayesian_optimization/af_optimizers/AFOptimizer.cpp
		optimizers/bayesian_optimization/af_optimizers/AFPSO.cpp
		#
		# algorithms: df-tr
		optimizers/trust_region/TrustRegionOptimization.cpp
		optimizers/trust_region/TrustRegionModel.cpp
		optimizers/trust_region/TrustRegionMath.cpp
		optimizers/ensemble_exp_value.cpp
		#
		# algorithms: hybrid
		hybrid_optimizer.cpp
		#
		# solvers
		../ThirdParty/snopt/handlers/SNOPTHandler.cpp
		../ThirdParty/snopt/handlers/SNOPTLoader.c
		../ThirdParty/snopt/handlers/LibraryHandler.c
		solvers/SNOPTSolver.cpp
		)

SET(OPTIMIZATION_TESTS
		# resources
		tests/test_resource_cases.h
		tests/test_resource_optimizer.h
		tests/test_resource_test_functions.h
		#
		# objective / case / optimizer
		# tests/objective/test_weightedsum.cpp
		# tests/test_case.cpp
		# tests/test_case_handler.cpp
		# tests/test_case_transfer_object.cpp
		# tests/test_normalizer.cpp
		# #
		# # constraints / constraint-handling
		# tests/constraints/test_bhp_constraint.cpp
		# tests/constraints/test_rate_constraint.cpp
		# tests/constraints/test_reservoir_boundary.cpp
		# tests/constraints/test_interwell_distance.cpp
		# tests/constraints/test_spline_well_length.cpp
		# tests/constraints/test_pseudo_cont_boundary_2d.cpp
		# tests/constraints/test_constraint_handler.cpp
		#
		# algorithms: pattern, stochastic
		# tests/optimizers/test_compass_search.cpp
		# tests/optimizers/test_apps.cpp
		# tests/optimizers/test_ga.cpp
		# tests/optimizers/test_pso.cpp
		# tests/optimizers/test_cma_es.cpp
		# tests/optimizers/test_vfsa.cpp
		# tests/optimizers/test_spsa.cpp
		# #
		# # algorithms: bayesian
		# tests/optimizers/test_ego.cpp
		#
		# algorithms: df-tr
		# tests/optimizers/test_tr-dfo.cpp
		tests/optimizers/test_en-tr-dfo.cpp
		tests/optimizers/test_tr-dfo_exp-value.cpp
		tests/optimizers/test_tr-model-data.hpp
		tests/optimizers/test_tr-support.hpp
		)
