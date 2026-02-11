#pragma once

#include "monomial.hpp"

namespace groebner_basis {

class LexicographicOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const {
        size_t max_size = std::max(lhs.getNumVariables(), rhs.getNumVariables());
        
        for (size_t i = 0; i < max_size; ++i) {
            auto lhs_degree = lhs.getDegreeOf(i);
            auto rhs_degree = rhs.getDegreeOf(i);
            
            if (lhs_degree != rhs_degree) {
                return lhs_degree > rhs_degree;
            }
        }
        
        return false;
    }
};

class GradedLexOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const {
        auto lhs_total = lhs.getTotalDegree();
        auto rhs_total = rhs.getTotalDegree();
        
        if (lhs_total != rhs_total) {
            return lhs_total > rhs_total;
        }
        
        return LexicographicOrder{}(lhs, rhs);
    }
};

class GradedReverseLexOrder {
public:
    bool operator()(const Monomial& lhs, const Monomial& rhs) const {
        auto lhs_total = lhs.getTotalDegree();
        auto rhs_total = rhs.getTotalDegree();
        
        if (lhs_total != rhs_total) {
            return lhs_total > rhs_total;
        }
        
        size_t max_size = std::max(lhs.getNumVariables(), rhs.getNumVariables());
        
        for (size_t i = max_size; i > 0; --i) {
            auto lhs_degree = lhs.getDegreeOf(i - 1);
            auto rhs_degree = rhs.getDegreeOf(i - 1);
            
            if (lhs_degree != rhs_degree) {
                return lhs_degree < rhs_degree;
            }
        }
        
        return false;
    }
};

}  // namespace groebner_basis
