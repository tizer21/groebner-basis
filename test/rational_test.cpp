#include <gtest/gtest.h>

#include "groebner_basis/fields/rational.hpp"

using namespace groebner_basis::fields;

class RationalTest : public ::testing::Test {
protected:
    void SetUp() override {
        zero = Rational(0);
        one = Rational(1);
        half = Rational(1, 2);
        third = Rational(1, 3);
        two_thirds = Rational(2, 3);
        negative_half = Rational(-1, 2);
        three_quarters = Rational(3, 4);
    }

    Rational zero;
    Rational one;
    Rational half;
    Rational third;
    Rational two_thirds;
    Rational negative_half;
    Rational three_quarters;
};

TEST_F(RationalTest, DefaultConstruction) {
    Rational r;
    EXPECT_EQ(r, zero);
}

TEST_F(RationalTest, ConstructionWithNumeratorOnly) {
    Rational r(5);
    EXPECT_EQ(r, Rational(5, 1));
    
    Rational r_neg(-3);
    EXPECT_EQ(r_neg, Rational(-3, 1));
}

TEST_F(RationalTest, ConstructionWithNumeratorAndDenominator) {
    Rational r(2, 3);
    EXPECT_EQ(r, two_thirds);
    
    Rational r2(4, 6);
    EXPECT_EQ(r2, two_thirds);
}

TEST_F(RationalTest, ConstructionReduction) {
    Rational r1(2, 4);
    EXPECT_EQ(r1, half);
    
    Rational r2(3, 9);
    EXPECT_EQ(r2, third);
    
    Rational r3(15, 25);
    EXPECT_EQ(r3, Rational(3, 5));
    
    Rational r4(-2, 4);
    EXPECT_EQ(r4, negative_half);
    
    Rational r5(2, -4);
    EXPECT_EQ(r5, negative_half);
    
    Rational r6(-2, -4);
    EXPECT_EQ(r6, half);
}

TEST_F(RationalTest, UnaryMinus) {
    EXPECT_EQ(-zero, zero);
    EXPECT_EQ(-one, Rational(-1));
    EXPECT_EQ(-half, negative_half);
    EXPECT_EQ(-negative_half, half);
}

TEST_F(RationalTest, Addition) {
    EXPECT_EQ(half + half, one);
    EXPECT_EQ(third + two_thirds, one);
    EXPECT_EQ(half + third, Rational(5, 6));
    
    EXPECT_EQ(zero + half, half);
    EXPECT_EQ(half + zero, half);
    
    EXPECT_EQ(half + negative_half, zero);
    EXPECT_EQ(negative_half + negative_half, Rational(-1));
}

TEST_F(RationalTest, AdditionAssignment) {
    Rational r = half;
    r += third;
    EXPECT_EQ(r, Rational(5, 6));
    
    r += zero;
    EXPECT_EQ(r, Rational(5, 6));
    
    r += negative_half;
    EXPECT_EQ(r, Rational(1, 3));
}

TEST_F(RationalTest, Subtraction) {
    EXPECT_EQ(one - half, half);
    EXPECT_EQ(two_thirds - third, third);
    EXPECT_EQ(half - third, Rational(1, 6));
    
    EXPECT_EQ(half - zero, half);
    EXPECT_EQ(zero - half, negative_half);
    
    EXPECT_EQ(half - negative_half, one);
}

TEST_F(RationalTest, SubtractionAssignment) {
    Rational r = one;
    r -= half;
    EXPECT_EQ(r, half);
    
    r -= third;
    EXPECT_EQ(r, Rational(1, 6));
    
    r -= zero;
    EXPECT_EQ(r, Rational(1, 6));
}

TEST_F(RationalTest, Multiplication) {
    EXPECT_EQ(half * half, Rational(1, 4));
    EXPECT_EQ(third * two_thirds, Rational(2, 9));
    EXPECT_EQ(half * two_thirds, Rational(1, 3));
    
    EXPECT_EQ(zero * half, zero);
    EXPECT_EQ(half * zero, zero);
    
    EXPECT_EQ(one * half, half);
    EXPECT_EQ(half * one, half);
    
    EXPECT_EQ(half * negative_half, Rational(-1, 4));
    EXPECT_EQ(negative_half * negative_half, Rational(1, 4));
}

TEST_F(RationalTest, MultiplicationAssignment) {
    Rational r = half;
    r *= third;
    EXPECT_EQ(r, Rational(1, 6));
    
    r *= one;
    EXPECT_EQ(r, Rational(1, 6));
    
    r *= zero;
    EXPECT_EQ(r, zero);
}

TEST_F(RationalTest, Division) {
    EXPECT_EQ(one / half, Rational(2));
    EXPECT_EQ(half / half, one);
    EXPECT_EQ(one / third, Rational(3));
    EXPECT_EQ(two_thirds / third, Rational(2));
    
    EXPECT_EQ(half / one, half);
    EXPECT_EQ(one / one, one);
    
    EXPECT_EQ(one / negative_half, Rational(-2));
    EXPECT_EQ(negative_half / negative_half, one);
}

TEST_F(RationalTest, DivisionAssignment) {
    Rational r = one;
    r /= half;
    EXPECT_EQ(r, Rational(2));
    
    r /= Rational(2);
    EXPECT_EQ(r, one);
    
    r /= third;
    EXPECT_EQ(r, Rational(3));
}

TEST_F(RationalTest, Equality) {
    EXPECT_EQ(zero, Rational(0));
    EXPECT_EQ(one, Rational(1));
    EXPECT_EQ(half, Rational(1, 2));
    EXPECT_EQ(two_thirds, Rational(2, 3));
    EXPECT_EQ(two_thirds, Rational(4, 6));
    
    EXPECT_NE(zero, one);
    EXPECT_NE(half, third);
    EXPECT_NE(half, negative_half);
}

TEST_F(RationalTest, LessThan) {
    EXPECT_LT(zero, one);
    EXPECT_LT(negative_half, zero);
    EXPECT_LT(negative_half, half);
    EXPECT_LT(third, half);
    EXPECT_LT(half, two_thirds);
    
    EXPECT_FALSE(one < zero);
    EXPECT_FALSE(half < third);
}

TEST_F(RationalTest, LessThanOrEqual) {
    EXPECT_LE(zero, one);
    EXPECT_LE(half, half);
    EXPECT_LE(negative_half, zero);
    EXPECT_LE(third, two_thirds);
    
    EXPECT_FALSE(one <= zero);
    EXPECT_FALSE(two_thirds <= third);
}

TEST_F(RationalTest, GreaterThan) {
    EXPECT_GT(one, zero);
    EXPECT_GT(half, third);
    EXPECT_GT(two_thirds, half);
    EXPECT_GT(zero, negative_half);
    
    EXPECT_FALSE(zero > one);
    EXPECT_FALSE(third > half);
}

TEST_F(RationalTest, GreaterThanOrEqual) {
    EXPECT_GE(one, zero);
    EXPECT_GE(half, half);
    EXPECT_GE(two_thirds, third);
    EXPECT_GE(zero, negative_half);
    
    EXPECT_FALSE(zero >= one);
    EXPECT_FALSE(third >= two_thirds);
}

TEST_F(RationalTest, OutputStream) {
    std::stringstream ss;
    
    ss << zero;
    EXPECT_EQ(ss.str(), "0");
    ss.str("");
    
    ss << one;
    EXPECT_EQ(ss.str(), "1");
    ss.str("");
    
    ss << half;
    EXPECT_EQ(ss.str(), "1/2");
    ss.str("");
    
    ss << negative_half;
    EXPECT_EQ(ss.str(), "-1/2");
    ss.str("");
    
    ss << Rational(4, 2);
    EXPECT_EQ(ss.str(), "2");
    ss.str("");
    
    ss << Rational(-3, 1);
    EXPECT_EQ(ss.str(), "-3");
}

TEST_F(RationalTest, ChainOperations) {
    Rational result = (half + third) * Rational(6, 5);
    EXPECT_EQ(result, one);
    
    result = (one - half) / third;
    EXPECT_EQ(result, Rational(3, 2));
    
    result = (-half) * (-two_thirds) / third;
    EXPECT_EQ(result, one);
}

TEST_F(RationalTest, NegativeNumbers) {
    EXPECT_EQ(Rational(-1, 2), negative_half);
    EXPECT_EQ(Rational(1, -2), negative_half);
    EXPECT_EQ(Rational(-1, -2), half);
    
    EXPECT_EQ(negative_half + negative_half, Rational(-1));
    EXPECT_EQ(negative_half * negative_half, Rational(1, 4));
    EXPECT_EQ(negative_half / negative_half, one);
}

TEST_F(RationalTest, MixedOperations) {
    EXPECT_EQ(half + third * two_thirds, Rational(13, 18));
    EXPECT_EQ((half + third) * two_thirds, Rational(5, 9));
    EXPECT_EQ(half / third + two_thirds, Rational(13, 6));
}

TEST_F(RationalTest, ReductionAfterOperations) {
    Rational r1 = Rational(2, 4) + Rational(2, 4);
    EXPECT_EQ(r1, one);
    
    Rational r2 = Rational(1, 3) * Rational(3, 1);
    EXPECT_EQ(r2, one);
    
    Rational r3 = Rational(10, 15) * Rational(3, 2);
    EXPECT_EQ(r3, one);
}
