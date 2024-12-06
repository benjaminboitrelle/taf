#include <boost/test/unit_test.hpp>

#include <sstream>

#define private public

#include "code/include/MAnalysis.h"

BOOST_AUTO_TEST_SUITE(tests_CreateGlobalResultDir)

BOOST_AUTO_TEST_CASE(test_CreateGlobalResultDir_init_not_done)
{
    auto analyse{MimosaAnalysis()};
    analyse.fInitDone = false;
    std::string message{""};
    BOOST_CHECK_EQUAL(message, analyse.CreateGlobalResultDir());
}

BOOST_AUTO_TEST_CASE(test_CreateGlobalResultDir_MimosaType_not_set)
{
    auto analyse{MimosaAnalysis()};
    analyse.fInitDone = true;
    std::string message{""};
    BOOST_CHECK_EQUAL(message, analyse.CreateGlobalResultDir());
}

BOOST_AUTO_TEST_CASE(test_CreateGlobalResultDir_pass)
{
    auto analyse{MimosaAnalysis()};
    analyse.fInitDone = true;
    BOOST_CHECK_EQUAL(true, analyse.CheckIfDone("init"));
    analyse.MimosaType = 4;
    analyse.fWorkingDirectory = "/Users/ben/Projects/cpp/taf/tests";
    auto created_dir = analyse.CreateGlobalResultDir();
    BOOST_CHECK_EQUAL("/Users/ben/Projects/cpp/taf/tests/results_ana_M4/", created_dir);
}
BOOST_AUTO_TEST_SUITE_END()