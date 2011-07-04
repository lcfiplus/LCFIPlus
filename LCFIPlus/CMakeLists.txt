########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( LCFIPlus )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( ROOT 5.28 REQUIRED )
FIND_PACKAGE( LCFIVertex 0.5.1 REQUIRED )

FOREACH( pkg Marlin MarlinUtil ROOT LCFIVertex )
    IF( ${pkg}_FOUND )
	INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
	LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
	ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()


INCLUDE_DIRECTORIES( ${LCFIVertex_DIR}/vertex_lcfi ${LCFIVertex_DIR}/boost )
LINK_LIBRARIES( -lEve )




### LIBRARY #################################################################

# include directories
INCLUDE_DIRECTORIES( "./include" )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )
ADD_DEFINITIONS( "-Wno-strict-aliasing" ) # avoid warnings in dict.cc

# add library
AUX_SOURCE_DIRECTORY( src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

