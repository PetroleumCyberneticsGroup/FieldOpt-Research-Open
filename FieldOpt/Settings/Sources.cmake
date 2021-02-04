SET(SETTINGS_HEADERS
	ensemble.h
	global.h
	model.h
	optimizer.h
	paths.h
	settings.h
	settings_exceptions.h
	simulator.h
	trajectory_importer.h
	helpers.hpp
	)

SET(SETTINGS_SOURCES
	ensemble.cpp
	global.cpp
	model.cpp
	optimizer.cpp
	paths.cpp
	settings.cpp
	simulator.cpp
	trajectory_importer.cpp
	)

SET(SETTINGS_TESTS
	tests/test_resource_example_file_paths.hpp
	tests/test_resource_multispline_wells_settings.hpp
	tests/test_resource_settings.hpp
	tests/test_settings.cpp
	tests/test_settings_model.cpp
	tests/test_settings_optimizer.cpp
	tests/test_settings_simulator.cpp
	tests/test_trajectory_importer.cpp
	)