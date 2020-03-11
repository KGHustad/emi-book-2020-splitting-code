add_library(ionic_model_grandi_cpu SHARED grandi_cpu.c)

if(OpenMP_C_FOUND)
    target_link_libraries(ionic_model_grandi_cpu PUBLIC OpenMP::OpenMP_C)
endif()
