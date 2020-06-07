#include "BigUIntBase.h"

namespace big {

bool BigUIntBase::isWellFormed() const {
    if (not isCorrectlySized()) { return false; }
    for (const auto& it : m_digits) {
        if (it & s_highBits) { return false; }
    }
    return true;
}

void BigUIntBase::resizeToFitVector(std::vector<size_t>& digits) {
    auto it = digits.rbegin();
    for (; it != digits.rend() && *it == 0ul; ++it)
        ;

    digits.resize(static_cast<size_t>(std::distance(it, digits.rend())));
    if (digits.empty()) { digits = {0ul}; }
}

void BigUIntBase::copyViaIterators(rlcIterator begin, rlcIterator end, rlIterator resultIt) {
    for (; begin != end; ++begin, ++resultIt) { *resultIt = *begin; }
}

} // namespace big