add_library(linear_solver_viennacl_cpu SHARED solve_viennacl_cpu.cpp)
target_include_directories(linear_solver_viennacl_cpu PUBLIC ../external)


if(OpenMP_CXX_FOUND)
    target_compile_definitions(linear_solver_viennacl_cpu PUBLIC VIENNACL_WITH_OPENMP)
    target_link_libraries(linear_solver_viennacl_cpu PUBLIC OpenMP::OpenMP_CXX)
endif()
