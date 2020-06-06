#include "BigUIntBase.h"

namespace big {

    bool BigUIntBase::isCorrectlySized() const {
        assert(digitCount() > 0);
        if (mostSignificantDigit() == 0ul) {
            return digitCount() == 1; // Is the value zero
        }
        return true;
    }

    bool BigUIntBase::isWellFormed() const {
        if (not isCorrectlySized()) {
            return false;
        }
        for (const auto &it : m_digits) {
            if (it & s_highBits) {
                return false;
            }
        }
        return true;
    }

    void BigUIntBase::resizeToFitVector(std::vector<size_t> &digits) {
        auto it = digits.rbegin();
        for (; it != digits.rend() && *it == 0ul; ++it);

        digits.resize(static_cast<size_t>(std::distance(it, digits.rend())));
        if (digits.empty()) {
            digits = {0ul};
        }
    }

    std::vector<size_t> BigUIntBase::shiftedCopy(size_t shiftAmount) const {
        std::vector<size_t> result(digitCount() + shiftAmount);
        std::copy(m_digits.begin(), m_digits.end(), result.begin() + shiftAmount);
        return result;
    }

} // namespace big