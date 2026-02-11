#include <gtest/gtest.h>
#include <sstream>
#include <vector>

#include "groebner_basis/monomial.hpp"

using namespace groebner_basis;

TEST(MonomialTest, DefaultConstructor) {
    Monomial m;
    EXPECT_EQ(m.getNumVariables(), 0u);
    EXPECT_EQ(m.getTotalDegree(), 0u);
    EXPECT_TRUE(m.isOne());
    std::ostringstream oss;
    oss << m;
    EXPECT_EQ(oss.str(), "1");
}

TEST(MonomialTest, ConstructorFromInitializerList) {
    Monomial m({1, 0, 2});
    EXPECT_EQ(m.getNumVariables(), 3u);
    EXPECT_EQ(m.getDegreeOf(0), 1u);
    EXPECT_EQ(m.getDegreeOf(1), 0u);
    EXPECT_EQ(m.getDegreeOf(2), 2u);
    EXPECT_EQ(m.getDegreeOf(3), 0u);
    EXPECT_EQ(m.getTotalDegree(), 3u);
    EXPECT_FALSE(m.isOne());
}

TEST(MonomialTest, ConstructorFromVector) {
    std::vector<unsigned int> exponents = {1, 0, 2};
    Monomial m(exponents);
    EXPECT_EQ(m.getNumVariables(), 3u);
    EXPECT_EQ(m.getDegreeOf(0), 1u);
    EXPECT_EQ(m.getDegreeOf(1), 0u);
    EXPECT_EQ(m.getDegreeOf(2), 2u);
    EXPECT_EQ(m.getTotalDegree(), 3u);
    EXPECT_FALSE(m.isOne());
}

TEST(MonomialTest, TrailingZerosTrimmed) {
    Monomial m({1, 0, 0});
    EXPECT_EQ(m.getNumVariables(), 1u);
    EXPECT_EQ(m.getDegreeOf(0), 1u);
    EXPECT_EQ(m.getTotalDegree(), 1u);

    Monomial m2({0, 0, 0});
    EXPECT_EQ(m2.getNumVariables(), 0u);
    EXPECT_TRUE(m2.isOne());

    Monomial m3({0, 2, 0});
    EXPECT_EQ(m3.getNumVariables(), 2u);
    EXPECT_EQ(m3.getDegreeOf(0), 0u);
    EXPECT_EQ(m3.getDegreeOf(1), 2u);
}

TEST(MonomialTest, Multiplication) {
    Monomial m1({1, 0, 2});
    Monomial m2({0, 3, 1});
    Monomial product = m1 * m2;
    EXPECT_EQ(product.getNumVariables(), 3u);
    EXPECT_EQ(product.getDegreeOf(0), 1u);
    EXPECT_EQ(product.getDegreeOf(1), 3u);
    EXPECT_EQ(product.getDegreeOf(2), 3u);
    EXPECT_EQ(product.getTotalDegree(), 7u);

    Monomial m3({1});
    Monomial m4({0, 3, 1});
    Monomial prod2 = m3 * m4;
    EXPECT_EQ(prod2.getNumVariables(), 3u);
    EXPECT_EQ(prod2.getDegreeOf(0), 1u);
    EXPECT_EQ(prod2.getDegreeOf(1), 3u);
    EXPECT_EQ(prod2.getDegreeOf(2), 1u);
}

TEST(MonomialTest, MultiplicationInPlace) {
    Monomial m1({1, 0});
    Monomial m2({0, 3});
    m1 *= m2;
    EXPECT_EQ(m1.getNumVariables(), 2u);
    EXPECT_EQ(m1.getDegreeOf(0), 1u);
    EXPECT_EQ(m1.getDegreeOf(1), 3u);
}

TEST(MonomialTest, Division) {
    Monomial dividend({2, 1, 3});
    Monomial divisor({1, 0, 1});
    Monomial quotient = dividend / divisor;
    EXPECT_EQ(quotient.getNumVariables(), 3u);
    EXPECT_EQ(quotient.getDegreeOf(0), 1u);
    EXPECT_EQ(quotient.getDegreeOf(1), 1u);
    EXPECT_EQ(quotient.getDegreeOf(2), 2u);
    EXPECT_EQ(quotient.getTotalDegree(), 4u);

    Monomial div2({1, 3});
    Monomial divs2({0, 3});
    Monomial q2 = div2 / divs2;
    EXPECT_EQ(q2.getNumVariables(), 1u);
    EXPECT_EQ(q2.getDegreeOf(0), 1u);
}

TEST(MonomialTest, DivisionInPlace) {
    Monomial m({2, 1});
    Monomial d({1, 0});
    m /= d;
    EXPECT_EQ(m.getDegreeOf(0), 1u);
    EXPECT_EQ(m.getDegreeOf(1), 1u);
}

TEST(MonomialTest, Divides) {
    Monomial m1({1, 0});
    Monomial m2({1, 1});
    Monomial m3({0, 0});
    Monomial m4({2, 0});

    EXPECT_TRUE(m1.divides(m1));
    EXPECT_TRUE(m1.divides(m2));
    EXPECT_TRUE(m3.divides(m1));
    EXPECT_FALSE(m4.divides(m1));
    EXPECT_TRUE(m1.divides(m4));

    Monomial m5({1});
    Monomial m6({1, 2});
    EXPECT_TRUE(m5.divides(m6));
    Monomial m7({2});
    EXPECT_FALSE(m7.divides(m6));
}

TEST(MonomialTest, LeastCommonMultiple) {
    Monomial m1({1, 0, 2});
    Monomial m2({0, 3, 1});
    Monomial lcm12 = ComputeLeastCommonMultiple(m1, m2);
    EXPECT_EQ(lcm12.getNumVariables(), 3u);
    EXPECT_EQ(lcm12.getDegreeOf(0), 1u);
    EXPECT_EQ(lcm12.getDegreeOf(1), 3u);
    EXPECT_EQ(lcm12.getDegreeOf(2), 2u);
    EXPECT_EQ(lcm12.getTotalDegree(), 6u);

    Monomial unit;
    Monomial lcm_unit = ComputeLeastCommonMultiple(m1, unit);
    EXPECT_EQ(lcm_unit.getDegreeOf(0), 1u);
    EXPECT_EQ(lcm_unit.getDegreeOf(1), 0u);
    EXPECT_EQ(lcm_unit.getDegreeOf(2), 2u);

    Monomial m3({1});
    Monomial m4({0, 3, 1});
    Monomial lcm34 = ComputeLeastCommonMultiple(m3, m4);
    EXPECT_EQ(lcm34.getNumVariables(), 3u);
    EXPECT_EQ(lcm34.getDegreeOf(0), 1u);
    EXPECT_EQ(lcm34.getDegreeOf(1), 3u);
    EXPECT_EQ(lcm34.getDegreeOf(2), 1u);
}

TEST(MonomialTest, Equality) {
    Monomial m1({1, 0, 2});
    Monomial m2({1, 0, 2});
    Monomial m3({1, 1, 2});
    Monomial m4({1, 0});
    EXPECT_TRUE(m1 == m2);
    EXPECT_FALSE(m1 == m3);
    EXPECT_FALSE(m1 == m4);
    EXPECT_TRUE(m1 != m3);
    EXPECT_FALSE(m1 != m2);

    Monomial empty;
    Monomial empty2;
    EXPECT_TRUE(empty == empty2);
}

TEST(MonomialTest, OutputStream) {
    Monomial unit;
    std::ostringstream oss1;
    oss1 << unit;
    EXPECT_EQ(oss1.str(), "1");

    Monomial m_single({1});
    std::ostringstream oss2;
    oss2 << m_single;
    EXPECT_EQ(oss2.str(), "x_0");

    Monomial m_power({0, 3});
    std::ostringstream oss3;
    oss3 << m_power;
    EXPECT_EQ(oss3.str(), "x_1^3");

    Monomial m_multi({1, 0, 2});
    std::ostringstream oss4;
    oss4 << m_multi;
    EXPECT_EQ(oss4.str(), "x_0*x_2^2");

    Monomial m_trim({0, 2, 0});
    std::ostringstream oss5;
    oss5 << m_trim;
    EXPECT_EQ(oss5.str(), "x_1^2");
}
