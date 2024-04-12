# CMake generated Testfile for 
# Source directory: /home/mhgk4/SCONE
# Build directory: /home/mhgk4/SCONE/Build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(Unit_Tests "/home/mhgk4/SCONE/Build/unitTests")
set_tests_properties(Unit_Tests PROPERTIES  _BACKTRACE_TRIPLES "/home/mhgk4/SCONE/CMakeLists.txt;251;add_test;/home/mhgk4/SCONE/CMakeLists.txt;0;")
add_test(Integration_Tests "/home/mhgk4/SCONE/Build/integrationTests")
set_tests_properties(Integration_Tests PROPERTIES  WORKING_DIRECTORY "/home/mhgk4/SCONE" _BACKTRACE_TRIPLES "/home/mhgk4/SCONE/CMakeLists.txt;252;add_test;/home/mhgk4/SCONE/CMakeLists.txt;0;")
subdirs("RandomNumbers")
subdirs("LinearAlgebra")
subdirs("SharedModules")
subdirs("Visualisation")
subdirs("ParticleObjects")
subdirs("NamedGrids")
subdirs("NuclearData")
subdirs("Geometry")
subdirs("Tallies")
subdirs("CollisionOperator")
subdirs("TransportOperator")
subdirs("UserInterface")
subdirs("PhysicsPackages")
subdirs("DataStructures")
