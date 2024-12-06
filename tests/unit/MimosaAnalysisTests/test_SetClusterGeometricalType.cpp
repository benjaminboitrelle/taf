#include <boost/test/unit_test.hpp>

#include <sstream>

#define private public

#include "code/include/MAnalysis.h"

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
