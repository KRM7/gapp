from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.build import check_min_cppstd
from conan.tools.files import copy

class GappRecipe(ConanFile):
    name = "gapp"
    version = "0.3.0"
    description = "A genetic algorithms library in C++ for single- and multi-objective optimization."
    topics = ("optimization", "multi-objective-optimization", "constrained-optimization", "genetic-algorithm", "cpp20")
    url = "https://github.com/KRM7/gapp"
    license = "MIT"

    settings = "os", "compiler", "build_type", "arch"

    options = { "shared": [True, False] }
    default_options = { "shared": False }

    exports_sources = "src/*", "CMakeLists.txt", "gapp.natvis"
    exports = "LICENSE"

    def layout(self):
        cmake_layout(self)

    def validate(self):
        check_min_cppstd(self, "20")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["GAPP_BUILD_TESTS"] = False
        tc.variables["GAPP_USE_LTO"] = False
        tc.variables["GAPP_USE_WERROR"] = False
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        copy(self, "LICENSE", src=self.recipe_folder, dst=self.package_folder)

    def package_info(self):
        self.cpp_info.libdirs = ["lib/" + str(self.settings.build_type)]
        self.cpp_info.libs = ["gapp"]
