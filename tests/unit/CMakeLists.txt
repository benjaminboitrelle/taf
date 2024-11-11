# Récupérer les fichiers sources de test
file(GLOB TEST_SOURCES "*.cpp")

# Pour chaque fichier de test, créer un exécutable
foreach(test_source ${TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} PRIVATE Boost::unit_test_framework)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()