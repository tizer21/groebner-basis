#pragma once

#include <cassert>
#include <ostream>

namespace groebner_basis::fields {
namespace detail {

inline bool IsPrime(long long number) {
    if (number == 1) {
        return false;
    }
    for (long long i = 2; i * i <= number; ++i) {
        if (number % i == 0) {
            return false;
        }
    }
    return true;
}

inline long long ComputeModularExponentiation(long long number, long long pow, long long modulus) {
    long long result = 1;
    while (pow > 0) {
        if (pow % 2 == 1) {
            result = (result * number) % modulus;
        }
        number = (number * number) % modulus;
        pow /= 2;
    }
    return result;

    if (pow == 0) {
        return 1;
    }
    if (pow % 2 == 0) {
        auto x = ComputeModularExponentiation(number, pow / 2, modulus);
        return x * x % modulus;
    }
    return ComputeModularExponentiation(number, pow - 1, modulus) * number % modulus;
}

}  // namespace detail

template <long long modulus>
class PrimeField {
    static_assert(modulus > 1, "Modulus must be greater than 1");
    static_assert(detail::IsPrime(modulus), "Modulus must be prime");

    using Integer = long long;

public:
    PrimeField() = default;

    explicit PrimeField(Integer value) : value_(reduce(value)) {}

    PrimeField operator-() const {
        return PrimeField(modulus - value_);
    }

    PrimeField& operator+=(const PrimeField& other) {
        value_ += other.value_;
        if (value_ >= modulus) {
            value_ -= modulus;
        }

        return *this;
    }

    PrimeField& operator-=(const PrimeField& other) {
        value_ -= other.value_;
        if (value_ < 0) {
            value_ += modulus;
        }

        return *this;
    }

    PrimeField& operator*=(const PrimeField& other) {
        value_ *= other.value_;
        value_ %= modulus;

        return *this;
    }

    PrimeField& operator/=(const PrimeField& other) {
        assert(other.value_ != 0 && "The divisor must not be zero");

        return *this *= other.getInversed();
    }

    friend PrimeField operator+(PrimeField lhs, const PrimeField& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend PrimeField operator-(PrimeField lhs, const PrimeField& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend PrimeField operator*(PrimeField lhs, const PrimeField& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend PrimeField operator/(PrimeField lhs, const PrimeField& rhs) {
        assert(rhs.value_ != 0 && "The divisor must not be zero");

        lhs /= rhs;
        return lhs;
    }

    friend bool operator==(const PrimeField& lhs, const PrimeField& rhs) {
        return lhs.value_ == rhs.value_;
    }

    friend bool operator!=(const PrimeField& lhs, const PrimeField& rhs) {
        return lhs.value_ != rhs.value_;
    }

    friend bool operator<(const PrimeField& lhs, const PrimeField& rhs) {
        return lhs.value_ < rhs.value_;
    }

    friend bool operator<=(const PrimeField& lhs, const PrimeField& rhs) {
        return lhs.value_ <= rhs.value_;
    }

    friend bool operator>(const PrimeField& lhs, const PrimeField& rhs) {
        return lhs.value_ > rhs.value_;
    }

    friend bool operator>=(const PrimeField& lhs, const PrimeField& rhs) {
        return lhs.value_ >= rhs.value_;
    }

    friend std::ostream& operator<<(std::ostream& out, const PrimeField& number) {
        out << number.value_;
        return out;
    }

private:
    static Integer reduce(Integer value) {
        value %= modulus;

        if (value < 0) {
            value += modulus;
        }

        return value;
    }

    PrimeField getInversed() const {
        return PrimeField(detail::ComputeModularExponentiation(value_, modulus - 2, modulus));
    }

    Integer value_ = 0;
};

}  // namespace groebner_basis::fields
