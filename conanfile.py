from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
from conan.tools.files import copy
import os


class SymPPConan(ConanFile):
    name = "sympp"
    version = "1.1.0"
    license = "BSD-3-Clause"
    author = "Leonardo Araujo dos Santos"
    url = "https://github.com/leonardoaraujosantos/SymPP"
    description = (
        "Modern C++20 symbolic mathematics library — a clean-room port of SymPy "
        "algorithms with SymPy wired in as the validation oracle."
    )
    topics = ("symbolic-math", "cas", "computer-algebra", "cpp20", "sympy")

    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": True, "fPIC": True}

    # GMP and MPFR are the hard runtime dependencies.
    requires = ("gmp/6.3.0", "mpfr/4.2.1")

    exports_sources = (
        "CMakeLists.txt", "cmake/*", "include/*", "src/*", "LICENSE",
    )

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        # Header/template-heavy C++ library: no C/C++ std settings to drop.
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.variables["SYMPP_BUILD_SHARED"] = self.options.shared
        tc.variables["SYMPP_BUILD_TESTS"] = False  # tests need Python+SymPy at configure
        tc.variables["SYMPP_BUILD_EXAMPLES"] = False
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        copy(self, "LICENSE", src=self.source_folder,
             dst=os.path.join(self.package_folder, "licenses"))

    def package_info(self):
        self.cpp_info.libs = ["sympp"]
        self.cpp_info.set_property("cmake_file_name", "SymPP")
        self.cpp_info.set_property("cmake_target_name", "SymPP::sympp")
        self.cpp_info.requires = ["gmp::gmp", "mpfr::mpfr"]
