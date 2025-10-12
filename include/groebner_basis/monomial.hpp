#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <numeric>
#include <ostream>
#include <vector>

namespace groebner_basis {

class Monomial {
    using Degree = unsigned int;
    using DegreeContainer = std::vector<Degree>;
    using Index = size_t;

public:
    Monomial() = default;

    explicit Monomial(std::initializer_list<Degree> degrees) : degrees_(degrees) {
        deleteTrailingZeros(&degrees_);
    }

    explicit Monomial(std::vector<Degree> degrees) : degrees_(std::move(degrees)) {
        deleteTrailingZeros(&degrees_);
    }

    bool divides(const Monomial& other) const {
        if (getNumVariables() > other.getNumVariables()) {
            return false;
        }
        for (Index idx = 0; idx != getNumVariables(); ++idx) {
            if (degrees_[idx] > other.degrees_[idx]) {
                return false;
            }
        }
        return true;
    }
    
    Degree getTotalDegree() const {
        return std::accumulate(degrees_.begin(), degrees_.end(), Degree(0));
    }

    Degree getDegreeOf(Index idx) const {
        return (idx < degrees_.size()) ? degrees_[idx] : Degree(0);
    }

    Index getNumVariables() const {
        return degrees_.size();
    }

    bool isOne() const {
        return degrees_.empty();
    }

    Monomial& operator*=(const Monomial& other) {
        Index other_size = other.degrees_.size();
        if (degrees_.size() < other_size) {
            degrees_.resize(other_size);
        }
        for (Index idx = 0; idx != other_size; ++idx) {
            degrees_[idx] += other.degrees_[idx];
        }

        return *this;
    }

    Monomial& operator/=(const Monomial& other) {
        assert(other.divides(*this) && "Monomial division requires that divisor divides dividend");

        for (Index idx = 0; idx != other.degrees_.size(); ++idx) {
            degrees_[idx] -= other.degrees_[idx];
        }
        deleteTrailingZeros(&degrees_);

        return *this;
    }

    friend bool operator==(const Monomial&, const Monomial&);
    friend bool operator!=(const Monomial&, const Monomial&);

    friend std::ostream& operator<<(std::ostream&, const Monomial&);

private:
    static void deleteTrailingZeros(std::vector<Degree>* degrees) {
        assert(degrees);

        auto new_size = degrees->size();
        while (new_size > 0 && (*degrees)[new_size - 1] == 0) {
            --new_size;
        }
        degrees->resize(new_size);
    }

    DegreeContainer degrees_;
};

inline Monomial operator*(Monomial lhs, const Monomial& rhs) {
    lhs *= rhs;
    return lhs;
}

inline Monomial operator/(Monomial lhs, const Monomial& rhs) {
    assert(rhs.divides(lhs) && "Monomial division requires that divisor divides dividend");

    lhs /= rhs;
    return lhs;
}

inline bool operator==(const Monomial& lhs, const Monomial& rhs) {
    return (lhs.degrees_ == rhs.degrees_);
}

inline bool operator!=(const Monomial& lhs, const Monomial& rhs)  {
    return !(lhs == rhs);
}

inline std::ostream& operator<<(std::ostream& os, const Monomial& monom) {
    if (monom.isOne()) {
        return os << '1';
    }

    bool is_first = true;
    for (size_t idx = 0; idx < monom.getNumVariables(); ++idx) {
        if (monom.degrees_[idx] == 0) {
            continue;
        }

        if (!is_first) {
            os << '*';
        }
        os << "x_" << idx;
        if (monom.degrees_[idx] > 1) {
            os << '^' << monom.degrees_[idx];
        }

        is_first = false;
    }

    return os;
}

inline Monomial ComputeLeastCommonMultiple(const Monomial& lhs, const Monomial& rhs) {
    auto num_variables = std::max(lhs.getNumVariables(), rhs.getNumVariables());
    std::vector<unsigned int> degrees(num_variables);
    for (size_t i = 0; i != num_variables; ++i) {
        degrees[i] = std::max(lhs.getDegreeOf(i), rhs.getDegreeOf(i));
    }

    return Monomial(degrees);
}

}  // namespace groebner_basis
