# Installing additional dependencies for LCFIPlus

## JSON parser by nlohmann
This JSON parser is used by the weaver interface code taken from LCCAnalyses.

```
git clone https://github.com/nlohmann/json.git nlohmann_json
```

## ONNX Runtime C/C++ library

The ONNX Runtime library is needed to perform the inference using the trained weights output from weaver.
The following instructions were adapted from [this link](https://stackoverflow.com/questions/63420533/setting-up-onnx-runtime-on-ubuntu-20-04-c-api).

### Download and install library
Here we assume local home directory installation.
```
mkdir /tmp/onnxInstall
cd /tmp/onnxInstall
wget -O onnx_archive.nupkg https://www.nuget.org/api/v2/package/Microsoft.ML.OnnxRuntime/1.19.2
unzip onnx_archive.nupkg
cp runtimes/linux-x64/native/libonnxruntime.so ~/.local/lib/
mkdir -p ~/.local/include/onnxruntime
cp -r build/native/include/*.h ~/.local/include/onnxruntime/
```

### Create CMake files

Additional cmake files are needed for cmake to find onnxruntime as a package.

```
mkdir -p ~/.local/share/cmake/onnxruntime
cd ~/.local/share/cmake/onnxruntime
cat <<'EOF' > onnxruntime-config.cmake
# Custom cmake config file by jcarius to enable find_package(onnxruntime) without modifying LIBRARY_PATH and LD_LIBRARY_PATH
#
# This will define the following variables:
#   onnxruntime_FOUND        -- True if the system has the onnxruntime library
#   onnxruntime_INCLUDE_DIRS -- The include directories for onnxruntime
#   onnxruntime_LIBRARIES    -- Libraries to link against
#   onnxruntime_CXX_FLAGS    -- Additional (required) compiler flags

include(FindPackageHandleStandardArgs)

# Assume we are in <install-prefix>/share/cmake/onnxruntime/onnxruntimeConfig.cmake
get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(onnxruntime_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

set(onnxruntime_INCLUDE_DIRS ${onnxruntime_INSTALL_PREFIX}/include)
set(onnxruntime_LIBRARIES onnxruntime)
set(onnxruntime_CXX_FLAGS "") # no flags needed


find_library(onnxruntime_LIBRARY onnxruntime
    PATHS "${onnxruntime_INSTALL_PREFIX}/lib"
)

add_library(onnxruntime SHARED IMPORTED)
set_property(TARGET onnxruntime PROPERTY IMPORTED_LOCATION "${onnxruntime_LIBRARY}")
set_property(TARGET onnxruntime PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${onnxruntime_INCLUDE_DIRS}")
set_property(TARGET onnxruntime PROPERTY INTERFACE_COMPILE_OPTIONS "${onnxruntime_CXX_FLAGS}")

find_package_handle_standard_args(onnxruntime DEFAULT_MSG onnxruntime_LIBRARY onnxruntime_INCLUDE_DIRS)
EOF
cat <<'EOF' > onnxruntime-version.cmake
# Custom cmake version file by jcarius

set(PACKAGE_VERSION "1.19.2")

# Check whether the requested PACKAGE_FIND_VERSION is compatible
if("${PACKAGE_VERSION}" VERSION_LESS "${PACKAGE_FIND_VERSION}")
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
  set(PACKAGE_VERSION_COMPATIBLE TRUE)
  if("${PACKAGE_VERSION}" VERSION_EQUAL "${PACKAGE_FIND_VERSION}")
    set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
```

## Add object file to MARLIN_DLL
```
export MARLIN_DLL=$MARLIN_DLL:$HOME/.local/lib/libonnxruntime.so
```
