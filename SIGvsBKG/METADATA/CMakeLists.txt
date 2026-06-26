# Set the minimum required CMake version:
cmake_minimum_required( VERSION 3.7 FATAL_ERROR )

# The name of the package:
project(UserAnalysis VERSION 1.0.0)


# Add ROOT system directory and require ROOT.
find_package( ROOT 6.28.00 REQUIRED COMPONENTS Core Hist RIO Tree ROOTDataFrame ROOTVecOps)


# Pick up a local version of the AtlasCMake code if it exists:
find_package( AtlasCMake QUIET )

# Find the project that we depend on:
find_package( AnalysisBase )

# Set up CTest:
atlas_ctest_setup()

# Set up a work directory project:
atlas_project( WorkDir 21.2.189
   USE AnalysisBase 21.2.189
    )

# Set up the runtime environment setup script(s):
lcg_generate_env( SH_FILE ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh )
install( FILES ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh
   DESTINATION . )


# Declare this as an ATLAS package
#atlas_subdir(MyMetaReader)

# Declare the package name:
atlas_subdir( MyMetadataReader )

# Component(s) in the package:
atlas_add_component( MyMetadataReader
                     src/*.cxx
                     LINK_LIBRARIES GaudiKernel AthenaBaseComps )




# Set up CPack:
atlas_cpack_setup()
