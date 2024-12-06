
set(TEST_SOURCES
    test_CheckIfDone.cpp
    test_SetClusterGeometricalType.cpp
)

add_executable(MimosaAnalysisTests)
target_sources(MimosaAnalysisTests PRIVATE ${TEST_SOURCES})
target_link_libraries(MimosaAnalysisTests PRIVATE Boost::unit_test_framework ${libname}) 
add_test(NAME MimosaAnalysisTests COMMAND MimosaAnalysisTests)