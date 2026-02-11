include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
    OVERRIDE_FIND_PACKAGE
)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(groebner_basis_tests
    monomial_test.cpp
    rational_test.cpp
)

target_link_libraries(groebner_basis_tests
    groebner-basis
    gtest_main
)

include(GoogleTest)
gtest_discover_tests(groebner_basis_tests)
