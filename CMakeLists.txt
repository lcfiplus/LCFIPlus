########################################################
# cmake file for building LCFIPlus package
# @author Tomohiko Tanabe (ICEPP, The University of Tokyo)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5 FATAL_ERROR)
########################################################


# project name
PROJECT( LCFIPlus )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 11 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

# Prevent CMake falls back to the latest standard the compiler does support
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# ONNX support option
OPTION( ENABLE_ONNX "Enable ONNX Runtime support for ML inference" OFF )

FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
#FIND_PACKAGE( ROOT REQUIRED COMPONENTS Minuit2 TMVA TreePlayer )
FIND_PACKAGE( ROOT REQUIRED COMPONENTS Minuit2 TMVA TreePlayer Gui Geom Eve Minuit XMLIO RGL Ged EG MLP )
#FIND_PACKAGE( Torch REQUIRED )
FIND_PACKAGE( LCFIVertex REQUIRED )

# Optional ONNX dependencies
IF( ENABLE_ONNX )
  MESSAGE( STATUS "ONNX support enabled" )
  FIND_PACKAGE( onnxruntime 1.17.1 REQUIRED)
  FIND_PACKAGE( nlohmann_json REQUIRED)
  ADD_DEFINITIONS( -DENABLE_ONNX )
ELSE()
  MESSAGE( STATUS "ONNX support disabled" )
ENDIF()

#set(JSON_BuildTests OFF CACHE INTERNAL "")
#add_subdirectory(nlohmann_json)

INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
ADD_DEFINITIONS ( ${Marlin_DEFINITIONS} )
INCLUDE_DIRECTORIES( SYSTEM ${MarlinUtil_INCLUDE_DIRS} )
ADD_DEFINITIONS ( ${MarlinUtil_DEFINITIONS} )
INCLUDE_DIRECTORIES( SYSTEM ${ROOT_INCLUDE_DIRS} )
INCLUDE_DIRECTORIES( SYSTEM ${LCFIVertex_INCLUDE_DIRS} )
#INCLUDE_DIRECTORIES( SYSTEM ${TORCH_INCLUDE_DIRS} )
# left here for backwards compatibility
INCLUDE_DIRECTORIES( SYSTEM ${LCFIVertex_ROOT}/vertex_lcfi ${LCFIVertex_ROOT}/boost )

SET( ROOT_DICT_CINT_DEFINITIONS -noIncludePaths -inlineInputHeader -I${PROJECT_SOURCE_DIR}/include )
#INCLUDE( ${ROOT_DICT_MACROS_FILE} )
INCLUDE( "${ILCUTIL_ROOT}/cmakemodules/MacroRootDict.cmake" )

SET( ROOT_DICT_INPUT_HEADERS
		${PROJECT_SOURCE_DIR}/include/lcfiplus.h
		${PROJECT_SOURCE_DIR}/include/EventStore.h
		${PROJECT_SOURCE_DIR}/include/LCIOStorer.h
		${PROJECT_SOURCE_DIR}/include/TreeStorer.h
		${PROJECT_SOURCE_DIR}/include/JetFinder.h
		${PROJECT_SOURCE_DIR}/include/process.h
		${PROJECT_SOURCE_DIR}/include/FlavorTag.h
		${PROJECT_SOURCE_DIR}/include/MakeNtuple.h
		${PROJECT_SOURCE_DIR}/include/TrainMVA.h
		${PROJECT_SOURCE_DIR}/include/ReadMVA.h
		${PROJECT_SOURCE_DIR}/include/testproc.h
		${PROJECT_SOURCE_DIR}/include/TrackNtuple.h
		${PROJECT_SOURCE_DIR}/include/VertexMassRecovery.h
		${PROJECT_SOURCE_DIR}/include/VertexNtuple.h
		${PROJECT_SOURCE_DIR}/include/vertexfinderdnn.h
		${PROJECT_SOURCE_DIR}/include/DNNProvider2.h
)

# Add ONNX-related headers if enabled
IF( ENABLE_ONNX )
  LIST( APPEND ROOT_DICT_INPUT_HEADERS
		${PROJECT_SOURCE_DIR}/include/MLInputGenerator.h
		${PROJECT_SOURCE_DIR}/include/MLMakeNtuple.h
		${PROJECT_SOURCE_DIR}/include/MLInferenceWeaver.h
  )
ENDIF()

LIST( APPEND ROOT_DICT_INPUT_HEADERS ${PROJECT_SOURCE_DIR}/include/LinkDef.h )
SET( ROOT_DICT_INCLUDE_DIRS ${Marlin_INCLUDE_DIRS} ${LCIO_INCLUDE_DIRS})
GEN_ROOT_DICT_SOURCES( dict.cc )
INSTALL(FILES ${ROOT_DICT_OUTPUT_DIR}/dict_rdict.pcm DESTINATION lib)
LIST(APPEND library_sources ${ROOT_DICT_OUTPUT_SOURCES} )

SET( binary_sources ./src/Driver.cc ./src/EventNavigator.cc )

### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )
	ADD_SUBDIRECTORY( ./doc )
ENDIF()

### LIBRARY #################################################################

# include directories
INCLUDE_DIRECTORIES( "./include" )
#LINK_DIRECTORIES( "./lib" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wno-effc++ -Wno-shadow" )
ADD_DEFINITIONS( "-Wno-long-long" )
ADD_DEFINITIONS( "-Wno-strict-aliasing" ) # avoid warnings in dict.cc

# shut up warnings in boost
#ADD_DEFINITIONS( "-Wno-unused-local-typedefs" )

# add library
SET( libname ${PROJECT_NAME} )
AUX_SOURCE_DIRECTORY( ./src library_sources )
LIST( REMOVE_ITEM library_sources ${binary_sources} )

# Remove ONNX-related source files if ONNX is disabled
IF( NOT ENABLE_ONNX )
  LIST( REMOVE_ITEM library_sources
    ./src/ONNXRuntime.cc
    ./src/MLInferenceWeaver.cc
    ./src/MLInputGenerator.cc
    ./src/MLMakeNtuple.cc
    ./src/WeaverInterface.cc
  )
ENDIF()

ADD_SHARED_LIBRARY( ${libname} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${libname} DESTINATION lib )

# Link libraries
SET( LINK_LIBS
	${Marlin_LIBRARIES}
	${ROOT_LIBRARIES}
	${ROOT_COMPONENT_LIBRARIES}
	${MarlinUtil_LIBRARIES}
	${LCFIVertex_LIBRARIES}
)

# Add ONNX libraries if enabled
IF( ENABLE_ONNX )
  LIST( APPEND LINK_LIBS
    nlohmann_json::nlohmann_json
    onnxruntime::onnxruntime
  )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname} ${LINK_LIBS} )

### EXECUTABLE MODE ########################################################

# make rootdict for executable
LIST( REMOVE_ITEM ROOT_DICT_INPUT_HEADERS ${PROJECT_SOURCE_DIR}/include/LinkDef.h )
LIST( APPEND ROOT_DICT_INPUT_HEADERS ${PROJECT_SOURCE_DIR}/include/Driver.h )
LIST( APPEND ROOT_DICT_INPUT_HEADERS ${PROJECT_SOURCE_DIR}/include/EventNavigator.h )
LIST( APPEND ROOT_DICT_INPUT_HEADERS ${PROJECT_SOURCE_DIR}/include/bin_LinkDef.h )
GEN_ROOT_DICT_SOURCES( bin_dict.cc )
LIST(APPEND binary_sources ${ROOT_DICT_OUTPUT_SOURCES} )

ADD_EXECUTABLE( lcfiplus_bin EXCLUDE_FROM_ALL ${binary_sources} )

SET_TARGET_PROPERTIES( lcfiplus_bin PROPERTIES COMPILE_FLAGS "-DBUILD_EVE" OUTPUT_NAME lcfiplus )
TARGET_LINK_LIBRARIES( lcfiplus_bin ${libname} torch ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

