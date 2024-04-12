# Install script for directory: /home/mhgk4/SCONE

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/mhgk4/SCONE/Build/RandomNumbers/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/LinearAlgebra/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/SharedModules/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/Visualisation/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/ParticleObjects/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/NamedGrids/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/NuclearData/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/Geometry/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/Tallies/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/CollisionOperator/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/TransportOperator/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/UserInterface/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/PhysicsPackages/cmake_install.cmake")
  include("/home/mhgk4/SCONE/Build/DataStructures/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/mhgk4/SCONE/Build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
