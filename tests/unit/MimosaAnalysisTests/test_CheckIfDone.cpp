#define BOOST_TEST_MODULE header - only multiunit test
// #include <boost/test/unit_test.hpp>

#include <boost/test/included/unit_test.hpp>
#include <sstream>

#define private public
#include "code/include/MAnalysis.h"

BOOST_AUTO_TEST_SUITE(test_CheckIfDone)

BOOST_AUTO_TEST_CASE(test_CheckIfDone_false_return)
{
    auto analyse{MimosaAnalysis()};

    analyse.fClearDone = false;

    auto check_options = {
        "clear",
        "init",
        "mimosall",
        "mimosapro",
        "mimosafakerate",
        "mimosacalibration",
        "mimosaminivectors",
        "mimosapro2planes",
        "mimosaimaging",
    };

    for (const auto &check : check_options)
    {
        BOOST_CHECK_EQUAL(false, analyse.CheckIfDone(check));
    }
}

BOOST_AUTO_TEST_CASE(test_CheckIfDone_multiple_inputs)
{
    auto analyse{MimosaAnalysis()};

    analyse.fInitDone = false;
    analyse.fMimosaProDone = false;

    auto check_options = "init,mimosapro";
    BOOST_CHECK_EQUAL(false, analyse.CheckIfDone(check_options));
}

BOOST_AUTO_TEST_SUITE_END()
