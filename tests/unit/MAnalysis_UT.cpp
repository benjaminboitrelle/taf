#define BOOST_TEST_MODULE MAnalysis_UT
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

BOOST_AUTO_TEST_SUITE(SetClusterGeometricalType)

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_OnePixel)
{
    MimosaAnalysis analysis;
    // Test case for one pixel
    analysis.SetClusterGeometricalType(1, 0, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType),
                      static_cast<int>(MAnalysis::ClusterType::one_pix));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_TwoPixelsRow)
{
    MimosaAnalysis analysis;

    // Test case for two pixels in a row
    analysis.SetClusterGeometricalType(2, 0, 2);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::two_pix_row));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_TwoPixelsColumn)
{
    MimosaAnalysis analysis;

    // Test case for two pixels in a column
    analysis.SetClusterGeometricalType(2, 2, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::two_pix_col));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_ThreePixelsL)
{
    MimosaAnalysis analysis;

    // Test case for three pixels forming an L-shape
    analysis.SetClusterGeometricalType(3, 2, 2);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::three_pix_L));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_ThreePixelsRow)
{
    MimosaAnalysis analysis;

    // Test case for three pixels in a row
    analysis.SetClusterGeometricalType(3, 3, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::three_pix_row));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_ThreePixelsColumn)
{
    MimosaAnalysis analysis;

    // Test case for three pixels in a column
    analysis.SetClusterGeometricalType(3, 0, 3);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::three_pix_col));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_FourPixelsSquare)
{
    MimosaAnalysis analysis;

    // Test case for four pixels in a square configuration
    analysis.SetClusterGeometricalType(4, 2, 2);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::four_pix_square));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_FourPixelsLRow)
{
    MimosaAnalysis analysis;

    // Test case for four pixels in a row L-shape
    analysis.SetClusterGeometricalType(4, 2, 3);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::four_pix_L_row));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_FourPixelsLCol)
{
    MimosaAnalysis analysis;

    // Test case for four pixels in a column L-shape
    analysis.SetClusterGeometricalType(4, 3, 2);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::four_pix_L_col));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_FourPixelsOthers)
{
    MimosaAnalysis analysis;

    // Test case for four pixels in any other configuration
    analysis.SetClusterGeometricalType(4, 0, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::four_others));
}

BOOST_AUTO_TEST_CASE(TestSetClusterGeometricalType_MoreThanFour)
{
    MimosaAnalysis analysis;

    // Test case for more than four pixels
    analysis.SetClusterGeometricalType(5, 3, 2);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::five_pix_squarerow));
    analysis.SetClusterGeometricalType(5, 2, 3);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::five_pix_squarecol));
    analysis.SetClusterGeometricalType(5, 0, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::five_others));
    analysis.SetClusterGeometricalType(6, 2, 3);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::six_pix_3col2row));
    analysis.SetClusterGeometricalType(6, 3, 2);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::six_pix_2col3row));
    analysis.SetClusterGeometricalType(6, 0, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::six_others));
    analysis.SetClusterGeometricalType(0, 0, 0);
    BOOST_CHECK_EQUAL(static_cast<int>(analysis.ClusterGeometricalType), static_cast<int>(MAnalysis::ClusterType::more_than_six));
}

BOOST_AUTO_TEST_SUITE_END()

// BOOST_AUTO_TEST_SUITE(GetFileNumber)
// BOOST_AUTO_TEST_CASE(TestGetFileNumber_noFile)
// {
//     auto analysis{MimosaAnalysis()};
//     const auto expected{0};
//     BOOST_CHECK_EQUAL(expected, analysis.GetFileNumber());
// }
// BOOST_AUTO_TEST_SUITE_END()

// Mocking ComputePixelPosition and ComputePixelPosition_UVToColRow for simplicity
// void ComputePixelPosition_UVToColRow(double u, double v, double &col, double &lin)
// {
//     // Simple mock: assuming u, v map directly to col, lin (you would adjust as needed for your logic)
//     col = u;
//     lin = v;
// }

// void ComputePixelPosition(int col, int lin, float &u, float &v)
// {
//     // Simple mock: assuming col, lin map directly back to u, v (you would adjust as needed)
//     u = static_cast<float>(col);
//     v = static_cast<float>(lin);
// }

BOOST_AUTO_TEST_SUITE(GetTrackDistantToClosestDiode)
BOOST_AUTO_TEST_CASE(Test_GetTrackDistantToClosestDiode)
{
    // Test case 1: Track exactly on a diode (distance should be 0)
    MimosaAnalysis analysis;
    constexpr auto initialDistance{1.0e+20f};
    SensorCoord3D<float> track_position{0.0f, 0.0f, 0.f}; // Assuming the track intersects exactly at (0, 0)

    // Assuming that the closest diode is at (0, 0)
    BOOST_CHECK_EQUAL(analysis.GetTrackDistantToClosestDiode(track_position), initialDistance);
}

BOOST_AUTO_TEST_SUITE_END()