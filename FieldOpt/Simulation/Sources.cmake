
SET(SIMULATION_HEADERS
	execution_scripts/execution_scripts.h
	results/adgprsresults.h
	results/eclresults.h
	results/results.h
	results/results_exceptions.h
	simulator_interfaces/adgprssimulator.h
	simulator_interfaces/driver_file_writers/adgprsdriverfilewriter.h
	simulator_interfaces/driver_file_writers/driver_parts/adgprs_driver_parts/adgprs_wellcontrols.h
	simulator_interfaces/driver_file_writers/driver_parts/adgprs_driver_parts/wellstre.h
	simulator_interfaces/driver_file_writers/driver_parts/driverpart.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/compdat.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/compsegs.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/ecldriverpart.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/schedule_section.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/wellcontrols.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welsegs.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welspecs.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/wsegvalv.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/schedule_insets.h
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/actionx.hpp
	simulator_interfaces/driver_file_writers/driver_parts/ix_driver_parts/flow_control_device.hpp
	simulator_interfaces/driver_file_writers/driver_parts/ix_driver_parts/ix_control.hpp
	simulator_interfaces/driver_file_writers/driver_parts/ix_driver_parts/report_tuning.hpp
	simulator_interfaces/driver_file_writers/ecldriverfilewriter.h
	simulator_interfaces/driver_file_writers/flowdriverfilewriter.h
	simulator_interfaces/driver_file_writers/ix_driver_file_writer.h
	simulator_interfaces/eclsimulator.h
	simulator_interfaces/flowsimulator.h
	simulator_interfaces/ix_simulator.h
	simulator_interfaces/simulator.h
	simulator_interfaces/simulator_exceptions.h
)

SET(SIMULATION_SOURCES
	results/adgprsresults.cpp
	results/eclresults.cpp
	simulator_interfaces/adgprssimulator.cpp
	simulator_interfaces/driver_file_writers/adgprsdriverfilewriter.cpp
	simulator_interfaces/driver_file_writers/driver_parts/adgprs_driver_parts/adgprs_wellcontrols.cpp
	simulator_interfaces/driver_file_writers/driver_parts/adgprs_driver_parts/wellstre.cpp
	simulator_interfaces/driver_file_writers/driver_parts/driverpart.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/compdat.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/compsegs.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/ecldriverpart.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/schedule_section.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/wellcontrols.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welsegs.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/welspecs.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/wsegvalv.cpp
	simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/schedule_insets.cpp
	simulator_interfaces/driver_file_writers/ecldriverfilewriter.cpp
	simulator_interfaces/driver_file_writers/flowdriverfilewriter.cpp
	simulator_interfaces/driver_file_writers/ix_driver_file_writer.cpp
	simulator_interfaces/eclsimulator.cpp
	simulator_interfaces/flowsimulator.cpp
	simulator_interfaces/ix_simulator.cpp
	simulator_interfaces/simulator.cpp
)

SET(SIMULATION_TESTS
	tests/test_resource_results.h
	tests/results/test_adgprsresults.cpp
	tests/results/test_eclresults.cpp
	tests/simulator_interfaces/driver_file_writers/adgprs_driver_file_writer.cpp
	tests/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/test_compdat.cpp
	tests/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/test_schedule_section.cpp
	tests/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/test_wellcontrols.cpp
	tests/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/test_welspecs.cpp
	tests/simulator_interfaces/driver_file_writers/driver_parts/ecl_driver_parts/test_schedule_inset.cpp
	tests/simulator_interfaces/driver_file_writers/flow_driver_file_writer.cpp
	tests/simulator_interfaces/test_adgprssimulator.cpp
	tests/simulator_interfaces/test_eclsimulator.cpp
	tests/simulator_interfaces/test_ix_simulator.cpp
)
