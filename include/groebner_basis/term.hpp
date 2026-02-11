#pragma once

#include <cassert>
#include <ostream>

#include "monomial.hpp"

namespace groebner_basis {

template <typename Coefficient>
class Term {
public:
    Term() = default;
    
    explicit Term(Coefficient coefficient) : coefficient_(coefficient) {}
    
    explicit Term(Monomial monomial) : monomial_(monomial) {}
    
    Term(Coefficient coefficient, Monomial monomial) 
        : coefficient_(coefficient), monomial_(monomial) {
        normalize(&coefficient_, &monomial_);
    }

    const Coefficient& getCoefficient() const {
        return coefficient_;
    }

    const Monomial& getMonomial() const {
        return monomial_;
    }

    bool isZero() const {
        return coefficient_ == Coefficient(0); 
    }

    bool divides(const Term& dividend) {
        if (isZero()) {
            return false;
        }
        if (dividend.isZero()) {
            return true;
        }
        return monomial_.divides(dividend.monomial_);
    }

    Term operator-() const {
        return Term(-coefficient_, monomial_);
    }

    Term& operator*=(const Term& other) {
        coefficient_ *= other.coefficient_;
        monomial_ *= other.monomial_;
        normalize(&coefficient_, &monomial_);

        return *this;
    }

    Term& operator/=(const Term& other) {
        assert(other.divides(*this) && "Term division requires that divisor divides dividend");

        if (!isZero()) {
            coefficient_ /= other.coefficient_;
            monomial_ /= other.monomial_;
        }

        return *this;
    }

    friend Term operator*(Term lhs, const Term& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend Term operator/(Term lhs, const Term& rhs) {
        assert(rhs.divides(lhs) && "Term division requires that divisor divides dividend");

        lhs /= rhs;
        return lhs;
    }

    bool operator==(const Term& other) const {
        return (coefficient_ == other.coefficient_) && (monomial_ == other.monomial_);
    }

    bool operator!=(const Term& other) const {
        return !(*this == other);
    }

    friend std::ostream& operator<<(std::ostream& os, const Term& term) {
        if (term.isZero() == 0) {
            return os << '0';
        }
        if (term.monomial_.isOne()) {
            return os << term.coefficient_;
        }

        return os << term.monomial_;
    }

private:
    static void normalize(Coefficient* coefficient, Monomial* monomial) {
        if (*coefficient == 0) {
            *monomial = Monomial();
        }
    }

    Coefficient coefficient_ = Coefficient(1);
    Monomial monomial_;
};

}  // namespace groebner_basis
