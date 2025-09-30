# Lib90B

Lib90B is a C++ library for estimating the entropy of data sources, derived from the [SP800-90B_EntropyAssessment](https://github.com/usnistgov/SP800-90B_EntropyAssessment) project by NIST. It implements statistical methods for assessing min-entropy as described in [NIST SP 800-90B](https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.800-90B.pdf). 

Lib90B currently, only offers the non-IID & IID test suites. Refactors of Restart and other tests are in progress.

## Features

- Non-IID min-entropy estimation suite  
- IID test suite
- Individual statistical tests: Most Common, Collision, Markov, Compression, LRS, t-Tuple, Multi-MCW, Lag, Multi-MMC, LZ78Y , Len_LRS, Chi Square, Permutation
- Support for 8-bit and arbitrary-width symbols  
- Optional unit tests via GoogleTest  


## License

NIST-developed software is provided as a public service. You may use, copy, and distribute copies of the software, including modified versions, provided that you retain this notice. The software is provided "AS IS" without warranty. Full license details are included in the source repository.


## Building Lib90B

Lib90B uses CMake and supports C++23. You can build it as a **static library**:

```bash
git clone https://github.com/HLRichardson-git/lib90b.git
cd lib90b
mkdir build && cd build
cmake -DLib90B_BUILD_TESTS=ON ..
cmake --build . --config Release
```

### Dependencies

libdivsufsort
 (for suffix array calculations, fetched automatically by CMake)

## Using Lib90B in your CMake project

Example CMakeLists.txt snippet to include Lib90B:

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)

include(FetchContent)

# Fetch Lib90B
FetchContent_Declare(
    lib90b
    GIT_REPOSITORY https://github.com/HLRichardson-git/lib90b.git
    GIT_TAG main
)
FetchContent_MakeAvailable(lib90b)

# Your executable
add_executable(my_app src/main.cpp)

# Link against Lib90B
target_link_libraries(my_app PRIVATE Lib90B)
target_include_directories(my_app PRIVATE ${lib90b_SOURCE_DIR}/include)
```

This approach will automatically download and make Lib90B available as part of your build, without requiring you to manually clone the repository into your project.

## Examples

Running the full Non-IID test suite

```cpp
#include <lib90b/non_iid.h>
#include <iostream>

int main() {
    lib90b::NonIidResult results = lib90b::nonIidTestSuite("path/to/data.bin");

    if (results.H_bitstring.has_value())
        std::cout << "H_bitstring: " << results.H_bitstring.value() << std::endl;
    if (results.H_original.has_value())
        std::cout << "H_original: " << results.H_original.value() << std::endl;
    if (results.min_entropy.has_value())
        std::cout << "Min entropy: " << results.min_entropy.value() << std::endl;

    return 0;
}
```

Running the full ID test suite

```cpp
#include <lib90b/iid.h>
#include <iostream>

int main() {
    lib90b::IidResult results = lib90b::iidTestSuite("path/to/data.bin");

    if (results.passed)
        std::cout << "IID Test Passed." << std::endl;

    return 0;
}
```

Running Individual Tests

```cpp
#include <lib90b/most_common.h>
#include <lib90b/t_tuple.h>

std::filesystem::path file_path = "entropy_samples/samples.bin";

// Most Common Test
auto most_common_result = lib90b::mostCommonTest(file_path);
std::cout << "Most common min entropy: " << most_common_result.h_original.value_or(0) << "\n";

// t-Tuple Test
auto ttuple_result = lib90b::tTupleTest(file_path, lib90b::SymbolMode::Bitstring);
std::cout << "t-Tuple min entropy: " << ttuple_result.h_bitstring.value_or(0) << "\n";
```

Note: Some tests can be called with `lib90b::SymbolMode::Bitstring` or `lib90b::SymbolMode::Original` to choose whether the test works on the raw symbols or their bitwise representation.
