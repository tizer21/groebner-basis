#pragma once

#include <cassert>
#include <numeric>
#include <ostream>

namespace groebner_basis::fields {

class Rational {
    using Integer = long long;

public:
    Rational() = default;

    explicit Rational(Integer numerator) : numerator_(numerator) {}
    
    Rational(Integer numerator, Integer denominator) : numerator_(numerator), denominator_(denominator) {
        assert(denominator_ != 0 && "The denominator must not be zero");
        reduce(&numerator_, &denominator_);
    }

    Rational operator-() const {
        return Rational(-numerator_, denominator_);
    }

    Rational& operator+=(const Rational& other) {
        numerator_ = (numerator_ * other.denominator_) +
                (other.numerator_ * denominator_);
        denominator_ *= other.denominator_;
        reduce(&numerator_, &denominator_);

        return *this;
    }

    Rational& operator-=(const Rational& other) {
        numerator_ = (numerator_ * other.denominator_) - (other.numerator_ * denominator_);
        denominator_ *= other.denominator_;
        reduce(&numerator_, &denominator_);

        return *this;
    }

    Rational& operator*=(const Rational& other) {
        numerator_ *= other.numerator_;
        denominator_ *= other.denominator_;
        reduce(&numerator_, &denominator_);

        return *this;
    }

    Rational& operator/=(const Rational& other) {
        assert(other != Rational(0) && "The divisor must not be zero");

        numerator_ *= other.denominator_;
        denominator_ *= other.numerator_;
        reduce(&numerator_, &denominator_);

        return *this;
    }

    friend bool operator==(const Rational&, const Rational&);
    friend bool operator!=(const Rational&, const Rational&);
    friend bool operator<(const Rational&, const Rational&);
    friend bool operator<=(const Rational&, const Rational&);
    friend bool operator>(const Rational&, const Rational&);
    friend bool operator>=(const Rational&, const Rational&);

    friend std::ostream& operator<<(std::ostream&, const Rational&);

private:
    static void reduce(Integer* numerator, Integer* denominator) {
        auto divisor = std::gcd(*numerator, *denominator);
        *numerator /= divisor;
        *denominator /= divisor;

        if (*denominator < 0) {
            *numerator = -(*numerator);
            *denominator = -(*denominator);
        }
    }

    Integer numerator_ = 0;
    Integer denominator_ = 1;
};

inline Rational operator+(Rational lhs, const Rational& rhs) {
    lhs += rhs;
    return lhs;
}

inline Rational operator-(Rational lhs, const Rational& rhs) {
    lhs -= rhs;
    return lhs;
}

inline Rational operator*(Rational lhs, const Rational& rhs) {
    lhs *= rhs;
    return lhs;
}

inline Rational operator/(Rational lhs, const Rational& rhs) {
    assert(rhs != Rational(0) && "The divisor must not be zero");

    lhs /= rhs;
    return lhs;
}

inline bool operator==(const Rational& lhs, const Rational& rhs) {
    return (lhs.numerator_ == rhs.numerator_) && (lhs.denominator_ == rhs.denominator_);
}

inline bool operator!=(const Rational& lhs, const Rational& rhs) {
    return !(lhs == rhs);
}

inline bool operator<(const Rational& lhs, const Rational& rhs) {
    return (lhs.numerator_ * rhs.denominator_) < (rhs.numerator_ * lhs.denominator_);
}

inline bool operator<=(const Rational& lhs, const Rational& rhs) {
    return (lhs < rhs) || (lhs == rhs);
}

inline bool operator>(const Rational& lhs, const Rational& rhs) {
    return !(lhs <= rhs);
}

inline bool operator>=(const Rational& lhs, const Rational& rhs) {
    return !(lhs < rhs);
}

inline std::ostream& operator<<(std::ostream &out, const Rational &rational) {
    out << rational.numerator_;
    if (rational.denominator_ != 1) {
        out << '/' << rational.denominator_;
    }
    return out;
}

} // namespace groebner_basis::fields
