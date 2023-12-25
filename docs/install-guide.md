# Install guide

[Requirements](#requirements)    
[Install with CMake](#install-with-cmake)  
[Install with vcpkg](#install-with-vcpkg)  
[Usage with CMake](#usage-with-cmake)  
[Using FetchContent](#using-fetchcontent)  
[CMake project options](#cmake-project-options)  
[Building the API docs](#building-the-api-docs)

## Requirements

You will need a C++20 compiler and CMake to build and use the library.  
There are no additional required dependencies, but [Catch2](https://github.com/catchorg/Catch2)
is needed to build the tests.

The full list of requirements are:

- C++20 compiler (gcc 11.0, clang 14.0, msvc 14.30 or later)
- CMake 3.21 or later
- Catch2 3.3 or later (optional, only needed for the tests)

The library works on Windows, Linux, and macOS. The tested compilers for
each of these platforms are:

- gcc and clang on Linux
- msvc and clang-cl on Windows
- gcc and clang (not AppleClang) on macOS

As the only real requirement for using the library is a compiler with C++20 support,
and the library doesn't include any platform or compiler specific code, other platforms
and compilers would probably also work, but only the ones listed above are tested.

Note that the standard library used also needs to support C++20 features, which means
that using libc\+\+ is not possible at this point. This is generally not a problem,
but on macOS the default standard library used by clang is libc\+\+ instead of libstdc\+\+,
which will cause build errors. In this case, the standard library should be manually
specified to be libstdc\+\+ instead.


## Install with CMake

The library uses CMake as its build system, which can also be used to install it:

```shell
# Clone the repository
git clone https://github.com/KRM7/gapp.git --branch v0.2.0
# Go to the library's build directory
cd gapp/build
# Configure cmake with the relevant options
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
# Build the library
cmake --build . --config Release
# Install the library
sudo cmake --install . --config Release
```

Alternatively, there is a utility script provided that can be used to
install the library in fewer steps:

```shell
git clone https://github.com/KRM7/gapp.git --branch v0.2.0
sudo bash gapp/build/install.sh
```

The script will perform the installation steps detailed above for all 3 configurations of the library.
If you want to configure CMake differently from the default, you can pass arguments to the
install script as you would do to CMake, for example:

```shell
sudo bash gapp/build/install.sh -DBUILD_TESTING=ON
```

### Install location

If you don't have elevated privileges, you will need to specify a location for the installation
using `CMAKE_INSTALL_PREFIX` during the configuration step:

```shell
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=/some/install/path
```

If you do this, you will also need to specify this path later on when using `find_package`:

```cmake
find_package(gapp REQUIRED PATH /some/install/path)
```

### Configurations

There are 3 configurations supported by the library: `Debug`, `Release`, and `RelWithDebInfo`.
If you don't specify the build type during the above installation steps, the `RelWithDebInfo` config
will be used as the default. If you want to install the library for more than 1 configuration,
you will need to repeat the above steps for all of the configurations you want.


### Running the tests

If Catch2 is installed, you can build and run the tests of the project as:

```shell
# In the build directory
cmake .. -DBUILD_TESTING=ON
cmake --build .
ctest --output-on-failure --schedule-random
```


## Install with vcpkg

The library can also be installed using Microsoft's
[vcpkg](https://github.com/microsoft/vcpkg/#getting-started) package manager:

```shell
vcpkg install gapp
```


## Usage with CMake

Once the library is installed, you can link against it using the namespaced CMake target
exported by the project. The name of this target is `gapp::gapp`.
A minimal CMakeLists.txt file of a project that depends on the library could for example
look like:

```cmake
project("example_project")

# Find the installed library
find_package(gapp CONFIG REQUIRED)

add_executable(example_project "example.cpp")
# Link the library to the executable
target_link_libraries(example_project PRIVATE gapp::gapp)
```

and the actual C++ code that uses the library could be:

```cpp
// example.cpp
#include <gapp/gapp.hpp> // include everything from the library

int main() {
    return gapp::rng::randomBool();
}
```


## Using FetchContent

The library can also be used as a subdirectory, allowing you to use CMake's `FetchContent`
instead of manually downloading and installing the library:

```cmake
project("example_project")

Include(FetchContent)

FetchContent_Declare(
    gapp
    GIT_REPOSITORY https://github.com/KRM7/gapp.git
    GIT_TAG        v0.2.0
)
FetchContent_MakeAvailable(gapp)

add_executable(example_project "example.cpp")
target_link_libraries(example_project PRIVATE gapp::gapp)
```

Note that if you use FetchContent instead of installing the library, the include paths of the library's
headers will be slightly different, ie. you will need to do `#include <gapp.hpp>` instead of the usual
`#include <gapp/gapp.hpp>`.

## CMake project options

- `BUILD_SHARED_LIBS` - When this option is `ON`, the library will be built as a shared library instead
    of a static one. Note that you have to choose between the static and shared versions of the library
    when installing it, you can't install both. The default value is `OFF`.

- `BUILD_TESTING` - When this option and the `GAPP_BUILD_TESTS` options are both `ON`, the tests will
    also be built along with the library. The default value is `ON`.

- `GAPP_BUILD_TESTS` - When this option and the `BUILD_TESTING` options are both `ON`, the tests will
    also be built along with the library. This option can be used to control when the tests are built
    when using the library as a subproject, without having to change the value of `BUILD_TESTING`.
    The default value is `ON`.

- `GAPP_BUILD_BENCHMARKS` - When this options and the `BUILD_TESTING` options are both `ON`, the
    benchmarks will also be built. The default value is `OFF`.

- `GAPP_BUILD_EXAMPLES` - When this option is `ON`, the usage examples in the examples folder will
    also be built. The default value is `OFF`.

- `GAPP_USE_WERROR` - When this option is `ON`, all warnings will be treated as errors during the
    build. The default value is `ON`.

- `GAPP_USE_LTO` - When this option is `ON`, the library will be built for using link-time optimization.
    The default value is `ON`.

- `GAPP_USE_MARCH_NATIVE` - When this option is `ON`, the library will be optimized specifically for the
    host architecture in the optimized configurations (Release and RelWithDebInfo). This has no effect on
    Debug builds. The default value is `OFF`.

- `GAPP_DISABLE_EXCEPTIONS` - When this option is `ON`, the library will be built with exception
    support disabled. The default value is `OFF`.

- `GAPP_DISABLE_RTTI` - When this option is `ON`, the library will be built without run-time type information.
    The default value is `OFF`.


The options are specified during the configuration step of CMake, for example:

```shell
# From the build directory
cmake .. -DGAPP_DISABLE_EXCEPTIONS=ON -DGAPP_DISABLE_RTTI=ON
cmake --build .
```

### Compiler options

If you want to specify additional compiler flags when building the library, you should do it through
the `GAPP_CXX_FLAGS` variable. The value of the `CMAKE_CXX_FLAGS` variable will be ignored by the project.
All of the flags specified through `GAPP_CXX_FLAGS` are in addition to the default flags,
they will not overwrite them.


## Building the API docs

The API documentation can be generated by running the `docs/api/generate_api_docs.sh` script.
The output will be in the `docs/api/sphinx-out` directory.
