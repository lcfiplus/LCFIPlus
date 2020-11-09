# Locates the tensorflow library and include dirs.

include(FindPackageHandleStandardArgs)
unset(Protobuf_FOUND)


find_path(Protobuf_INCLUDE_DIR
	NAMES
	google
	HINTS
	/home/goto/local/include
	/usr/include/tensorflow
	/usr/local/include/tensorflow)

# set Tensorflow_FOUND
find_package_handle_standard_args(Protobuf Protobuf_INCLUDE_DIR)

# set external variables for usage in CMakeLists.txt
if(Protobuf_FOUND)
	set(Protobuf_INCLUDE_DIRS ${Protobuf_INCLUDE_DIR})
endif()


# hide locals from GUI
mark_as_advanced(Protobuf_INCLUDE_DIRS)
