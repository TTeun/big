#include "BigUInt.h"

#include "BigInt.h"

#include <algorithm>
#include <bitset>
#include <iomanip>
#include <iostream>
#include <random>

namespace big {
static_assert(BigUInt::s_maxDigit <= std::numeric_limits<size_t>::max() / BigUInt::s_maxDigit,
              "s_base^2 should not exceed the maximum size_t");
static_assert(BigUInt::s_base % 2 == 0, "s_base should be even");

const size_t BigUIntBase::s_maxDigit; // So that it can be used in std::min in BigUInt::divisionSubRoutine

/***************** Some static functions *****************/
inline static size_t ceilingIntegerDivision(size_t a, size_t b) {
    return (a + b - 1) / b;
}

static size_t modPower(size_t base, size_t exponent, size_t mod) {
    // assert(base < BigUIntBase::s_base);
    if (exponent == 0ul) {
        return 1ul;
    } else if (exponent == 1ul) {
        return base % mod;
    }
    size_t aux(1ul);
    size_t copy(base);
    while (exponent > 1ul) {
        if (exponent % 2ul != 0ul) {
            aux *= copy;
            aux %= mod;
        }
        exponent /= 2ul;
        copy *= copy;
        copy %= mod;
    }
    copy *= aux;
    return copy % mod;
}

void bubbleViaIterators(rlIterator thisIt, const rlIterator& thisEnd) {
    uint32_t carry = 0ul;
    for (; thisIt != thisEnd; ++thisIt) {
        *thisIt += carry;
        if (*thisIt & BigUIntBase::s_highBits) {
            carry = static_cast<uint32_t>(BigUIntBase::divideByBase(*thisIt));
            *thisIt &= BigUIntBase::s_lowBits;
        } else {
            carry = 0ul;
        }
    }
}

/***************** Constructors *****************/
BigUInt::BigUInt(std::vector<size_t>&& digits, bool isAlreadyCorrectlySized) : BigUIntBase(std::move(digits)) {
    if (m_digits.empty()) { m_digits = {0ul}; }
    if (not isAlreadyCorrectlySized) { resizeToFit(); }
}

BigUInt::BigUInt(const std::string& val) {
    m_digits          = {0ul};
    size_t startIndex = val.length() % s_decimalsInDigit;
    size_t readBlock;
    if (startIndex != 0ul) {
        *this *= s_tenPowerInDigit;
        std::istringstream iss(val.substr(0, startIndex));
        iss >> readBlock;
        *this += readBlock;
    }
    while (startIndex + s_decimalsInDigit <= val.length()) {
        *this *= s_tenPowerInDigit;
        std::istringstream iss(val.substr(startIndex, s_decimalsInDigit));
        iss >> readBlock;
        *this += readBlock;
        startIndex += s_decimalsInDigit;
    }
}

/***************** Operators *****************/
/** Addition **/
BigUInt& BigUInt::operator+=(const BigUInt& rhs) {
    // ToDo make function to double numbers;
    if (this == &rhs) { return *this *= 2ul; }
    // assert(isWellFormed() && rhs.isWellFormed());
    if (rhs.isZero()) {
        return *this;
    } else if (isZero()) {
        return *this = rhs;
    } else if (rhs.digitCount() <= 2ul) {
        return *this += rhs.value();
    } else {
        resize(std::max(digitCount(), rhs.digitCount()) + 1ul);
        addViaIterators(rlBegin(), rhs.rlcBegin(), rhs.rlcEnd());
        reduceSizeByOneIfNeeded();
    }
    // assert(isWellFormed());
    return *this;
}

BigUInt& BigUInt::operator+=(size_t rhs) {
    if (rhs != 0ul) {
        if (isZero()) {
            return *this = rhs;
        } else if (rhs <= s_additionRoom) {
            resize(digitCount() + 1ul);
            carryAdditionViaIterators(rlBegin(), rhs);
            reduceSizeByOneIfNeeded();
        } else {
            *this += BigUInt(rhs);
        }
    }
    return *this;
}

BigUInt BigUInt::operator+(size_t rhs) const {
    return BigUInt(*this) += rhs;
}

BigUInt BigUInt::operator+(const BigUInt& rhs) const {
    auto copy = *this;
    return copy += rhs;
}

/** Subtraction **/
BigUInt& BigUInt::operator-=(const BigUInt& rhs) {
    subtractViaIterators(rlBegin(), rhs.rlcBegin(), rhs.rlcEnd());
    resizeToFit();
    return *this;
}

BigUInt BigUInt::operator-(const BigUInt& rhs) const {
    auto copy = *this;
    return copy -= rhs;
}

/** Multiplication **/
BigUInt& BigUInt::operator*=(const size_t rhs) {
    if (rhs == 0 || isZero()) {
        m_digits = {0ul};
    } else if (rhs != 1ul) {
        if (rhs < s_base) {
            multiplyBySingleDigit(rhs);
        } else {
            reserve(digitCount() + 2ul);
            *this *= BigUInt({rhs & s_lowBits, (rhs / s_base) & s_lowBits}, true);
        }
    }
    return *this;
}

BigUInt& BigUInt::operator*=(const BigUInt& rhs) {
    // assert(isWellFormed());
    // assert(rhs.isWellFormed());
    if (rhs.isZero()) {
        m_digits = {0ul};
    } else {
        if (this == &rhs) {
            // ToDo make function
            return *this *= BigUInt(rhs);
        }

        switch (rhs.digitCount()) {
            case 1ul: multiplyBySingleDigit(rhs.leastSignificantDigit()); break;
            default: *this = rhs.digitCount() > digitCount() ? multiply(*this, rhs) : multiply(rhs, *this); break;
        }
    }
    // assert(isWellFormed());
    return *this;
}

BigUInt BigUInt::operator*(size_t rhs) const {
    auto copy = *this;
    return copy *= rhs;
}

BigUInt BigUInt::operator*(const BigUInt& rhs) const {
    if (*this < rhs) {
        return multiply(*this, rhs);
    } else {
        return multiply(rhs, *this);
    }
}

/** Division **/
BigUInt& BigUInt::operator/=(size_t divisor) {
    if (divisor < s_base) {
        divideByLessThanBase(divisor);
    } else {
        *this = *this / divisor;
    }
    return *this;
}

BigUInt& BigUInt::operator/=(const BigUInt& divisor) {
    *this = *this / divisor;
    return *this;
}

BigUInt BigUInt::operator/(const BigUInt& divisor) const {
    // assert(divisor != 0ul);
    if (&divisor == this) { return BigUInt(1); }
    if (divisor == 1ul) { return BigUInt(*this); }
    if (divisor.digitCount() == 1ul && divisor.mostSignificantDigit() < s_base) {
        auto copy = *this;
        copy.divideByLessThanBase(divisor.mostSignificantDigit());
        return copy;
    }
    if (*this < divisor) { return 0; }
    const size_t m = this->digitCount();
    const size_t n = divisor.digitCount();

    if (m < n) { return 0ul; }
    if (m == 1) { return mostSignificantDigit() / divisor.mostSignificantDigit(); }
    BigUInt copy(*this);
    return longDivision(copy, divisor, false);
}

BigUInt BigUInt::operator/(size_t divisor) const {
    if (not(divisor & s_highBits)) {
        auto copy = *this;
        copy.divideByLessThanBase(divisor);
        return copy;
    } else {
        return *this / BigUInt(divisor);
    }
}

/** Modulo **/
BigUInt& BigUInt::operator%=(size_t mod) {
    *this = *this % mod;
    return *this;
}

BigUInt& BigUInt::operator%=(const BigUInt& mod) {
    // assert(mod != 0ul);
    if (*this < mod) {
        return *this;
    } else if (mod.digitCount() == 1ul) {
        return *this %= mod.mostSignificantDigit();
    } else {
        longDivision(*this, mod, true);
    }
    resizeToFit();
    return *this;
}

size_t BigUInt::operator%(size_t mod) const {
    if (mod == 2ul) { return leastSignificantDigit() % 2ul; }
    if (mod >= std::numeric_limits<size_t>::max() / mod) {
        auto copy = *this % BigUInt(mod);
        return value();
    }
    size_t result = 0ul;
    for (size_t i = 0; i != digitCount(); ++i) {
        result += (digitAt(i) % mod) * modPower(s_base % mod, i, mod);
        result %= mod;
    }
    return result;
}

BigUInt BigUInt::operator%(const BigUInt& mod) const {
    auto copy = *this;
    switch (mod.digitCount()) {
        case 1ul: return copy %= mod.leastSignificantDigit();
        case 2ul: return copy %= mod.mostSignificantDigit() * s_base + mod.leastSignificantDigit();
        default: return copy %= mod;
    }
}

/** Comparison **/
bool BigUInt::operator<=(const BigUInt& rhs) const {
    // assert(isWellFormed() && rhs.isCorrectlySized());

    if (digitCount() != rhs.digitCount()) { return digitCount() < rhs.digitCount(); }

    auto thisIt = lrcBegin();
    auto rhsIt  = rhs.lrcBegin();
    for (; thisIt != lrcEnd(); ++thisIt, ++rhsIt) {
        if (*thisIt != *rhsIt) { return *thisIt <= *rhsIt; }
    }
    return true;
}

/** Friends **/
BigUInt power(const BigUInt& base, size_t exponent) {
    switch (exponent) {
        case 0ul: return BigUInt(1);
        case 1ul: return BigUInt(base);
        default:
            BigUInt aux(1);
            BigUInt copy(base);

            while (exponent > 1ul) {
                if (exponent % 2ul == 1ul) { aux *= copy; }
                exponent /= 2ul;
                copy.square();
            }
            copy *= aux;

            return copy;
    }
}

/***************** Builders *****************/
BigUInt BigUInt::createRandom(size_t numberOfDigits) {
    // assert(numberOfDigits > 0ul);
    std::random_device                    rd;
    std::mt19937                          gen(rd());
    std::uniform_int_distribution<size_t> dis(0, s_maxDigit);

    BigUInt result;
    result.resize(numberOfDigits);
    for (size_t i = 0; i != numberOfDigits; ++i) { result.m_digits[i] = dis(gen); }
    result.resizeToFit();
    return result;
}

/***************** Output *****************/
std::string BigUInt::toString() const {
    std::stringstream ss;
    ss << *lrcBegin();
    for (auto it = lrcBegin() + 1ul; it != lrcEnd(); ++it) { ss << *it; }
    return ss.str();
}

std::string BigUInt::toDecimalString() const {
    // Terrible, don't use
    std::stringstream ss;
    auto              copy = *this;
    while (copy != 0ul) {
        size_t decimalDigit = 0ul;
        for (size_t i = 0ul; i != copy.digitCount(); ++i) { decimalDigit += (modPower(s_base % 10, i, 10ul) * copy.digitAt(i)) % 10ul; }
        decimalDigit %= 10;
        ss << std::setfill('0') << std::setw(1) << decimalDigit;
        copy /= 10;
    }
    auto reversedResult = ss.str();
    std::reverse(reversedResult.begin(), reversedResult.end());
    return reversedResult;
}

std::string BigUInt::toBinaryString() const {
    std::stringstream ss;
    ss << std::bitset<s_bitsPerDigit>(mostSignificantDigit());
    if (digitCount() > 1ul) {
        for (auto lrIt = lrcBegin() + 1ul; lrIt != lrcEnd(); ++lrIt) {
            ss << std::setfill('0') << std::setw(32) << std::bitset<s_bitsPerDigit>(*lrIt);
        }
    }
    auto str = ss.str();
    str.erase(0, std::min(str.find_first_not_of('0'), str.size() - 1));
    return str;
}

std::ostream& operator<<(std::ostream& os, const BigUInt& bigUnsignedInt) {
    os << "( ";
    for (auto it = bigUnsignedInt.lrcBegin(); it != bigUnsignedInt.lrcEnd(); ++it) { os << *it << " "; }
    os << ")_" << BigUInt::s_base;
    return os;
}

/***************** Internal *****************/
void BigUInt::bubble(size_t startIndex) {
    // assert(not m_digits.empty());
    // assert(startIndex < digitCount());
    // assert(isWellFormed());

    auto     it    = rlBegin() + startIndex;
    uint32_t carry = 0ul;
    for (; it != rlEnd(); ++it) {
        *it += carry;
        if (*it & s_highBits) {
            carry = static_cast<uint32_t>(divideByBase(*it));
            *it &= s_lowBits;
        } else {
            carry = 0ul;
        }
    }
    // assert(carry == 0ul);
    // assert(isWellFormed());
}

void BigUInt::divideByLessThanBase(size_t factor) {
    // assert(not(factor & s_highBits));
    uint32_t carry  = 0ul;
    auto     thisIt = lrBegin();
    for (; thisIt != lrEnd(); ++thisIt) {
        *thisIt += carry * s_base;
        carry = static_cast<uint32_t>(*thisIt % factor);
        *thisIt /= factor;
    }
    resizeToFit();
}

void BigUInt::square() {
    const auto copy = *this;
    *this *= copy;
    bubble(0);
}

void BigUInt::multiplyBySingleDigit(const size_t digit) {
    if (digitCount() == 1ul) {
        size_t val = digit * leastSignificantDigit();
        m_digits   = val & s_highBits ? std::vector{val & s_lowBits, divideByBase(val)} : std::vector{val};
    } else {
        resize(digitCount() + 1ul);
        multiplyBySingleDigitViaIterators(rlBegin(), rlEnd(), digit);
        reduceSizeByOneIfNeeded();
    }
}

/***************** Static helpers *****************/
/** Addition **/
void BigUInt::carryAdditionViaIterators(rlIterator resultIt, size_t carry) {
    // assert(carry != 0ul);
    if (carry == 1ul) {
        carryUnitAdditionViaIterators(resultIt);
    } else {
        // assert(carry + s_base < std::numeric_limits<size_t>::max());
        for (;; ++resultIt) {
            // assert(resultIt != resultEnd);
            *resultIt += carry;
            if (not(*resultIt & s_highBits)) { return; }
            carry = divideByBase(*resultIt);
            *resultIt &= s_lowBits;
        }
    }
}

void BigUInt::carryUnitAdditionViaIterators(rlIterator resultIt) {
    for (;; ++resultIt) {
        if (*resultIt < s_maxDigit) {
            ++*resultIt;
            return;
        }
        *resultIt = 0ul;
    }
}

void BigUInt::addViaIterators(rlIterator resultIt, rlcIterator rhsIt, rlcIterator rhsEnd) {
    static const uint32_t addBlockSize = 201u;

    uint8_t carry = 0u;
    size_t  tempValue;
    for (; rhsEnd - rhsIt >= addBlockSize;) {
        tempValue = *resultIt + *rhsIt + carry;
        *resultIt = tempValue & s_lowBits;
        ++resultIt;
        ++rhsIt;
        for (uint32_t dummy = 1u; dummy != addBlockSize; ++dummy) {
            tempValue = *resultIt + *rhsIt + (tempValue & s_highBits ? 1u : 0u);
            *resultIt = tempValue & s_lowBits;
            ++resultIt;
            ++rhsIt;
        }
        carry = (tempValue & s_highBits ? 1u : 0u);
    }

    switch (rhsEnd - rhsIt) {
        case 0: break;
        case 1:
            tempValue = *resultIt + carry + *rhsIt;
            *resultIt = tempValue & s_lowBits;
            ++resultIt;
            carry = (tempValue & s_highBits ? 1u : 0u);
            break;
        default:
            tempValue = *resultIt + carry + *rhsIt;
            *resultIt = tempValue & s_lowBits;
            ++resultIt;
            ++rhsIt;
            const auto limit = static_cast<uint32_t>(rhsEnd - rhsIt);
            for (uint32_t dummy = 0ul; dummy < limit; ++dummy) {
                tempValue = *resultIt + *rhsIt + (tempValue & s_highBits ? 1u : 0u);
                *resultIt = tempValue & s_lowBits;
                ++resultIt;
                ++rhsIt;
            }
            carry = (tempValue & s_highBits ? 1u : 0u);
            break;
    }
    if (carry) {
        // assert(carry == 1ul);
        carryUnitAdditionViaIterators(resultIt);
    }
}

void BigUInt::addMultipleViaIterators(rlIterator resultIt, rlcIterator rhsIt, rlcIterator rhsEnd, size_t multiplier) {
    // assert(multiplier != 0ul);
    static const uint32_t addBlockSize = 5u;

    uint32_t carry = 0ul;
    size_t   tempValue;
    for (; rhsEnd - rhsIt >= addBlockSize;) {
        tempValue = *resultIt + *rhsIt * multiplier + carry;
        *resultIt = tempValue & s_lowBits;
        ++resultIt;
        ++rhsIt;
        for (uint32_t dummy = 1u; dummy != addBlockSize; ++dummy) {
            tempValue = *resultIt + *rhsIt * multiplier + divideByBase(tempValue);
            *resultIt = tempValue & s_lowBits;
            ++resultIt;
            ++rhsIt;
        }
        carry = divideByBase(tempValue);
    }

    switch (rhsEnd - rhsIt) {
        case 0: break;
        case 1:
            tempValue = *resultIt + *rhsIt * multiplier + carry;
            *resultIt = tempValue & s_lowBits;
            ++resultIt;
            carry = divideByBase(tempValue);
            break;
        default:
            tempValue = *resultIt + *rhsIt * multiplier + carry;
            *resultIt = tempValue & s_lowBits;
            ++resultIt;
            ++rhsIt;
            const auto limit = static_cast<uint32_t>(rhsEnd - rhsIt);
            for (uint32_t dummy = 0ul; dummy < limit; ++dummy) {
                tempValue = *resultIt + *rhsIt * multiplier + divideByBase(tempValue);
                *resultIt = tempValue & s_lowBits;
                ++resultIt;
                ++rhsIt;
            }
            carry = divideByBase(tempValue);
            break;
    }
    switch (carry) {
        case 0: return;
        case 1: carryUnitAdditionViaIterators(resultIt); break;
        default: carryAdditionViaIterators(resultIt, carry); break;
    }
}

void BigUInt::subtractViaIterators(rlIterator thisIt, rlcIterator rhsIt, rlcIterator rhsEnd) {
    // assert(std::distance(thisIt, thisEnd) >= std::distance(rhsIt, rhsEnd));
    uint8_t carry = 0u;
    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        if (*thisIt >= *rhsIt + carry) {
            *thisIt -= *rhsIt + carry;
            carry = 0u;
        } else {
            *thisIt -= (*rhsIt + carry) - s_base;
            carry = 1u;
        }
    }
    if (carry != 0u) {
        while (*thisIt == 0ul) {
            // assert(thisIt != thisEnd);
            *thisIt = s_maxDigit;
            ++thisIt;
        }
        --*thisIt;
    }
}

/** Multiplication **/
BigUInt BigUInt::multiply(const BigUInt& smaller, const BigUInt& larger) {
    //    assert(smaller.digitCount() <= larger.digitCount());
    BigUInt result;
    result.resize(smaller.digitCount() + larger.digitCount());
    multiplySortedViaIterators(result.rlBegin(), smaller.rlcBegin(), smaller.rlcEnd(), larger.rlcBegin(), larger.rlcEnd());
    result.reduceSizeByTwoIfNeeded();
    return result;
}

void BigUInt::multiplyBySingleDigitViaIterators(rlIterator resultIt, const rlIterator resultEnd, const size_t rhs) {
    uint32_t carry = 0ul;
    for (; resultIt != resultEnd; ++resultIt) {
        *resultIt = *resultIt * rhs + carry;
        if (*resultIt & s_highBits) {
            carry = static_cast<uint32_t>(divideByBase(*resultIt));
            *resultIt &= s_lowBits;
        } else {
            carry = 0ul;
        }
    }
}

void BigUInt::karatsubaMultiplyViaIterators(
    rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd) {
    // assert((largeEnd - largeIt) >= (smallEnd - smallIt));
    const auto smallSize = static_cast<size_t>(smallEnd - smallIt);
    // assert(smallSize >= s_karatsubaLowerLimit);
    const auto   largeSize  = static_cast<size_t>(largeEnd - largeIt);
    const size_t splitIndex = smallSize / 2ul;

    BigUInt high1(largeIt + splitIndex, largeEnd);
    BigUInt high2(smallIt + splitIndex, smallIt + smallSize);

    BigUInt z0;
    z0.resize(2ul * splitIndex + 1ul);
    multiplySortedViaIterators(z0.rlBegin(), smallIt, smallIt + splitIndex, largeIt, largeIt + splitIndex);
    z0.resizeToFit();

    BigUInt z2;
    z2.resize(smallSize % 2ul + largeSize + 1ul);
    multiplySortedViaIterators(z2.rlBegin(), smallIt + splitIndex, smallEnd, largeIt + splitIndex, largeEnd);
    z2.resizeToFit();

    high1.resize(high1.digitCount() + 1ul);
    high2.resize(high2.digitCount() + 1ul);

    addViaIterators(high1.rlBegin(), largeIt, largeIt + splitIndex);
    addViaIterators(high2.rlBegin(), smallIt, smallIt + splitIndex);

    BigUInt z1;
    z1.resize(smallSize % 2ul + largeSize + 2ul);

    high1.resizeToFit();
    high2.resizeToFit();

    multiplyViaIterators(z1.rlBegin(), high1.rlcBegin(), high1.rlcEnd(), high2.rlcBegin(), high2.rlcEnd());

    subtractViaIterators(z1.rlBegin(), z2.rlcBegin(), z2.rlcEnd());
    subtractViaIterators(z1.rlBegin(), z0.rlcBegin(), z0.rlcEnd());

    copyViaIterators(z0.rlcBegin(), z0.rlcEnd(), resultIt);
    addViaIterators(resultIt + splitIndex, z1.rlcBegin(), z1.rlcEnd());
    addViaIterators(resultIt + 2ul * splitIndex, z2.rlcBegin(), z2.rlcEnd());
}

void BigUInt::splitOneMultiplicationViaIterators(
    rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd) {
    const auto largeSize = static_cast<size_t>(largeEnd - largeIt);
    const auto smallSize = static_cast<size_t>(smallEnd - smallIt);
    // assert(largeSize >= s_karatsubaLowerLimit);
    // assert(largeSize >= smallSize);
    const size_t splitIndex = largeSize / 2ul;

    multiplyViaIterators(resultIt, smallIt, smallEnd, largeIt, largeIt + splitIndex);
    BigUInt high;
    high.resize(largeSize - splitIndex + smallSize + 1ul);
    multiplyViaIterators(high.rlBegin(), largeIt + splitIndex, largeEnd, smallIt, smallEnd);
    high.resizeToFit();
    addViaIterators(resultIt + splitIndex, high.rlcBegin(), high.rlcEnd());
}

void BigUInt::multiplySortedViaIterators(
    rlIterator resultIt, rlcIterator smallIt, const rlcIterator smallEnd, rlcIterator largeIt, const rlcIterator largeEnd) {
    const auto smallSize = static_cast<size_t>(smallEnd - smallIt);
    const auto largeSize = static_cast<size_t>(largeEnd - largeIt);
    // assert(largeSize >= smallSize);

    if (largeSize < s_karatsubaLowerLimit) {
        schoolMultiply(resultIt, smallIt, smallEnd, largeIt, largeSize);
    } else if (smallSize < s_karatsubaLowerLimit) {
        splitOneMultiplicationViaIterators(resultIt, smallIt, smallEnd, largeIt, largeEnd);
    } else if (2ul * smallSize <= largeSize) {
        splitOneMultiplicationViaIterators(resultIt, smallIt, smallEnd, largeIt, largeEnd);
    } else if (smallSize < s_toomCook3LowerLimit) {
        karatsubaMultiplyViaIterators(resultIt, smallIt, smallEnd, largeIt, largeEnd);
    } else if (smallSize < s_toomCook4LowerLimit) {
        toomCook_3(resultIt, smallIt, smallEnd, largeIt, largeEnd);
    } else {
        toomCook_4(resultIt, smallIt, smallEnd, largeIt, largeEnd);
    }
}

void BigUInt::toomCook_3(rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd) {
    // assert((largeEnd - largeIt) >= (smallEnd - smallIt));
    const size_t i = (smallEnd - smallIt) / 3ul;

    const auto [m0, m1, m2] = static_cast<std::tuple<BigInt, BigInt, BigInt>>(splitThree(smallIt, smallEnd, i));
    const auto [n0, n1, n2] = static_cast<std::tuple<BigInt, BigInt, BigInt>>(splitThree(largeIt, largeEnd, i));

    const BigInt p_aux      = m0 + m2;
    const BigInt p_minusOne = p_aux - m1;
    const BigInt p_minusTwo = 2ul * (p_minusOne + m2) - m0;

    const BigInt q_aux      = n0 + n2;
    const BigInt q_minusOne = q_aux - n1;
    const BigInt q_minusTwo = 2ul * (q_minusOne + n2) - n0;

    const BigInt a0         = m0 * n0;
    const BigInt r_one      = (p_aux + m1) * (q_aux + n1);
    const BigInt r_minusOne = p_minusOne * q_minusOne;
    const BigInt r_minusTwo = q_minusTwo * p_minusTwo;
    const BigInt a4         = m2 * n2;

    BigInt a1;
    BigInt a2;
    BigInt a3;

    a3 = (r_minusTwo - r_one) / 3ul;
    a1 = (r_one - r_minusOne) / 2ul;
    a2 = r_minusOne - a0;
    a3 = (a2 - a3) / 2ul + 2ul * a4;
    a2 = a2 + a1 - a4;
    a1 -= a3;

    addViaIterators(resultIt, a0.magnitude().rlcBegin(), a0.magnitude().rlcEnd());
    addViaIterators(resultIt + i, a1.magnitude().rlcBegin(), a1.magnitude().rlcEnd());
    addViaIterators(resultIt + 2ul * i, a2.magnitude().rlcBegin(), a2.magnitude().rlcEnd());
    addViaIterators(resultIt + 3ul * i, a3.magnitude().rlcBegin(), a3.magnitude().rlcEnd());
    addViaIterators(resultIt + 4ul * i, a4.magnitude().rlcBegin(), a4.magnitude().rlcEnd());
}

void BigUInt::toomCook_4(rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd) {
    const size_t i              = (smallEnd - smallIt) / 4ul;
    const auto [m0, m1, m2, m3] = static_cast<std::tuple<BigInt, BigInt, BigInt, BigInt>>(splitFour(smallIt, smallEnd, i));
    auto [n0, n1, n2, n3]       = static_cast<std::tuple<BigInt, const BigInt, const BigInt, BigInt>>(splitFour(largeIt, largeEnd, i));

    BigInt n02 = n0 + n2;
    BigInt n13 = n1 + n3;
    BigInt m02 = m0 + m2;
    BigInt m13 = m1 + m3;

    const big::BigInt r_one      = (n02 + n13) * (m02 + m13);
    const big::BigInt r_minusOne = (n02 - n13) * (m02 - m13);
    n02 += 3ul * n2;
    n13 += n1 + 7ul * n3;
    m02 += 3ul * m2;
    m13 += m1 + 7ul * m3;
    const big::BigInt r_two      = (n02 + n13) * (m02 + m13);
    const big::BigInt r_minusTwo = (n02 - n13) * (m02 - m13);
    const big::BigInt r_three    = (n0 + 3ul * n1 + 9ul * n2 + 27ul * n3) * (m0 + 3ul * m1 + 9ul * m2 + 27ul * m3);

    n0 *= m0;
    n3 *= m3;

    BigInt a1 = (-20ll * n0 - 720ll * n3 + 60ll * r_one - 30ll * r_minusOne - 15ll * r_two + 3ll * r_minusTwo + 2ll * r_three) / 60ul;
    BigInt a2 = (-30ll * n0 + 96ll * n3 + 16ll * r_one + 16ll * r_minusOne - r_two - r_minusTwo) / 24ul;
    BigInt a3 = (10ll * n0 + 360ll * n3 - 14ll * r_one - r_minusOne + 7ll * r_two - r_minusTwo - r_three) / 24ul;
    BigInt a4 = (6ll * n0 - 120ll * n3 - 4ll * r_one - 4ll * r_minusOne + r_two + r_minusTwo) / 24ul;
    BigInt a5 = (-10ll * n0 - 360ll * n3 + 10ll * r_one + 5ll * r_minusOne - 5ll * r_two - r_minusTwo + r_three) / 120ul;

    addViaIterators(resultIt, n0.magnitude().rlcBegin(), n0.magnitude().rlcEnd());
    addViaIterators(resultIt + i, a1.magnitude().rlcBegin(), a1.magnitude().rlcEnd());
    addViaIterators(resultIt + 2ul * i, a2.magnitude().rlcBegin(), a2.magnitude().rlcEnd());
    addViaIterators(resultIt + 3ul * i, a3.magnitude().rlcBegin(), a3.magnitude().rlcEnd());
    addViaIterators(resultIt + 4ul * i, a4.magnitude().rlcBegin(), a4.magnitude().rlcEnd());
    addViaIterators(resultIt + 5ul * i, a5.magnitude().rlcBegin(), a5.magnitude().rlcEnd());
    addViaIterators(resultIt + 6ul * i, n3.magnitude().rlcBegin(), n3.magnitude().rlcEnd());
}

void BigUInt::schoolMultiply(rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, size_t largeSize) {
    for (size_t i = 0; i != largeSize; ++i) {
        if (*largeIt != 0ul) { addMultipleViaIterators(resultIt + i, smallIt, smallEnd, *largeIt); }
        ++largeIt;
    }
}

/** Division **/
size_t BigUInt::divisionSubRoutine(const lrcIterator leftToRightConstIt,
                                   const lrcIterator leftToRightConstEnd,
                                   const rlIterator  rightToLeftIt,
                                   const BigUInt&    divisor) {
    if (lessThanViaIterators(leftToRightConstIt, leftToRightConstEnd, divisor.lrcBegin(), divisor.lrcEnd())) { return 0ul; }
    // assert(divisor != 0ul);
    // assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);
    const size_t divisorSize  = divisor.digitCount();
    const auto   dividendSize = static_cast<size_t>(leftToRightConstEnd - leftToRightConstIt);
    // assert(dividendSize <= divisorSize + 1ul);
    // assert(dividendSize >= divisorSize);

    size_t correction = 0ul;
    while (not lessThanShiftedRhsViaIterators(leftToRightConstIt, leftToRightConstEnd, divisor.lrcBegin(), divisor.lrcEnd(), 1)) {
        subtractViaIterators(rightToLeftIt + 1ul, divisor.rlcBegin(), divisor.rlcEnd());
        correction += BigUInt::s_base;
    }

    size_t quotientEstimate;
    if (dividendSize == divisorSize) {
        quotientEstimate = (*leftToRightConstIt) / divisor.mostSignificantDigit();
    } else {
        quotientEstimate = (*leftToRightConstIt * BigUInt::s_base + *(leftToRightConstIt + 1ul)) / divisor.mostSignificantDigit();
    }
    quotientEstimate    = std::min(quotientEstimate, s_maxDigit);
    const size_t offset = *leftToRightConstIt == 0 ? 1ul : 0ul;

    BigUInt closestMultipleEstimate = createWithRoom(divisor.digitCount() + 1ul);
    closestMultipleEstimate         = divisor;
    closestMultipleEstimate.resize(closestMultipleEstimate.digitCount() + 1ul);
    multiplyBySingleDigitViaIterators(closestMultipleEstimate.rlBegin(), closestMultipleEstimate.rlEnd(), quotientEstimate);
    closestMultipleEstimate.resizeToFit();
    while (greaterThanViaIterators(
        closestMultipleEstimate.lrcBegin(), closestMultipleEstimate.lrcEnd(), leftToRightConstIt + offset, leftToRightConstEnd)) {
        --quotientEstimate;
        closestMultipleEstimate -= divisor;
    }
    subtractViaIterators(rightToLeftIt, closestMultipleEstimate.rlcBegin(), closestMultipleEstimate.rlcEnd());
    return quotientEstimate + correction;
}

BigUInt BigUInt::longDivision(BigUInt& u, BigUInt v, bool findRemainder) {
    // Knuth "Art of computer programming" volume 2
    const size_t d     = s_maxDigit / v.mostSignificantDigit();
    const size_t uSize = u.digitCount();
    const size_t n     = v.digitCount();
    const size_t m     = u.digitCount() - n;
    u.resize(uSize + 1ul);

    multiplyBySingleDigitViaIterators(u.rlBegin(), u.rlEnd(), d);
    multiplyBySingleDigitViaIterators(v.rlBegin(), v.rlEnd(), d);
    std::vector<size_t> quotient(m + 1ul, 0ul);
    for (long long j = m; j >= 0l; --j) {
        const size_t leadingU = (u.digitAt(n + j) << s_bitsPerDigit) + u.digitAt(n + j - 1);
        size_t       q        = leadingU / v.mostSignificantDigit();
        size_t       r        = leadingU % v.mostSignificantDigit();

        if ((q & s_base) || (q * v.digitAt(n - 2ul) > s_base * r + u.digitAt(n + j - 2))) {
            --q;
            r += v.mostSignificantDigit();
            if (r < s_base) {
                if ((q & s_base) || (q * v.digitAt(n - 2ul) > s_base * r + u.digitAt(n + j - 2))) {
                    --q;
                    r += v.mostSignificantDigit();
                }
            }
        }
        BigUInt temp = q * v;
        temp.resize(n + 1ul);
        if (lessThanViaIterators(u.lrcBegin() - m + j, u.lrcEnd() - j, temp.lrcBegin(), temp.lrcEnd())) {
            --q;
            temp -= v;
        }
        quotient[j] = q;
        subtractViaIterators(u.rlBegin() + j, temp.rlcBegin(), temp.rlcEnd());
    }
    if (findRemainder) {
        u.resize(n);
        u.divideByLessThanBase(d);
    }
    return BigUInt(std::move(quotient), false);
}

BigUInt BigUInt::longDivisionAfterAdjustingDivisor(BigUInt& dividend, const BigUInt& divisor) {
    // assert(divisor <= dividend);
    // assert(divisor.mostSignificantDigit() * 2ul >= BigUInt::s_base);
    assert(false);

    size_t       m = dividend.digitCount();
    const size_t n = divisor.digitCount();
    if (m <= n + 1) { return BigUInt::divisionSubRoutine(dividend.lrcBegin(), dividend.lrcEnd(), dividend.rlBegin(), divisor); }

    std::vector<size_t> divisorDigits(m - n + 1ul, 0ul);
    while (m > n + 1) {
        const size_t splitIndex = m - n - 1ul;
        divisorDigits[splitIndex] =
            BigUInt::divisionSubRoutine(dividend.lrcBegin(), dividend.lrcBegin() + n + 1, dividend.rlBegin() + splitIndex, divisor);
        dividend.resizeToFit();
        m = dividend.digitCount();
    }
    if (m > n || dividend.mostSignificantDigit() >= divisor.mostSignificantDigit()) {
        divisorDigits[0ul] = BigUInt::divisionSubRoutine(dividend.lrcBegin(), dividend.lrcEnd(), dividend.rlBegin(), divisor);
    }
    bubbleViaIterators(divisorDigits.begin(), divisorDigits.end());
    return BigUInt(std::move(divisorDigits), false);
}

/** Comparison **/
bool BigUInt::lessThanShiftedRhsViaIterators(
    lrcIterator thisIt, lrcIterator thisEnd, lrcIterator rhsIt, lrcIterator rhsEnd, size_t trailingZeroesOfRhs) {
    if (static_cast<size_t>(thisEnd - thisIt) != static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs) {
        return static_cast<size_t>(thisEnd - thisIt) < static_cast<size_t>(rhsEnd - rhsIt) + trailingZeroesOfRhs;
    }

    for (; rhsIt != rhsEnd; ++thisIt, ++rhsIt) {
        if (*thisIt != *rhsIt) { return *thisIt < *rhsIt; }
    }
    return false;
}

} // namespace big