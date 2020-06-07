#ifndef __BIG_U_INT__H__
#define __BIG_U_INT__H__

#include "BigUIntBase.h"

#include <cmath>
#include <ostream>
#include <tuple>

namespace big {

class BigUInt : public BigUIntBase {
public:
    static const size_t s_karatsubaLowerLimit = 200ul;
    static const size_t s_toomCook3LowerLimit = 800ul;
    static const size_t s_toomCook4LowerLimit = 2000ul;

    static BigUInt longDivision(BigUInt& dividend, const BigUInt& divisor);

public:
    /***************** Constructors *****************/
    BigUInt() { init(0); }
    BigUInt(size_t val) { init(val); }
    BigUInt(BigUInt&& other) noexcept : BigUIntBase(std::move(other.m_digits)) { }
    BigUInt(const BigUInt& other) : BigUIntBase(std::vector<size_t>(other.m_digits.begin(), other.m_digits.end())) {
        // assert(isWellFormed());
    }
    BigUInt(std::vector<size_t>&& digits, bool isAlreadyCorrectlySized);
    explicit BigUInt(const std::string& val);
    BigUInt(rlcIterator it, rlcIterator endIt) : BigUIntBase({it, endIt}) {
        // assert(isWellFormed());
    }

    /***************** Operators *****************/
    /** Assignment **/
    inline BigUInt& operator=(const BigUInt& rhs) {
        m_digits = rhs.m_digits;
        return *this;
    }
    inline BigUInt& operator=(size_t rhs) {
        init(rhs);
        return *this;
    }

    /** Addition **/
    BigUInt& operator+=(size_t rhs);
    BigUInt& operator+=(const BigUInt& rhs);
    BigUInt  operator+(size_t rhs) const;
    BigUInt  operator+(const BigUInt& rhs) const;

    /** Subtraction **/
    BigUInt& operator-=(const BigUInt& rhs);
    BigUInt  operator-(const BigUInt& rhs) const;

    /** Multiplication **/
    BigUInt& operator*=(size_t rhs);
    BigUInt& operator*=(const BigUInt& rhs);
    BigUInt  operator*(size_t rhs) const;
    BigUInt  operator*(const BigUInt& rhs) const;

    /** Division **/
    BigUInt& operator/=(size_t divisor);
    BigUInt& operator/=(const BigUInt& divisor);
    BigUInt  operator/(const BigUInt& divisor) const;
    BigUInt  operator/(size_t divisor) const;

    /** Modulo **/
    BigUInt& operator%=(size_t mod);
    BigUInt& operator%=(const BigUInt& mod);
    size_t   operator%(size_t mod) const;
    BigUInt  operator%(const BigUInt& mod) const;

    /** Comparison **/
    inline bool operator==(const BigUInt& rhs) const {
        //        assert(isWellFormed());
        //        assert(rhs.isWellFormed());
        return m_digits == rhs.m_digits;
    }
    bool operator!=(const BigUInt& rhs) const { return !(rhs == *this); }
    bool operator<(const BigUInt& rhs) const {
        //        assert(isWellFormed());
        //        assert(rhs.isWellFormed());
        return lessThanViaIterators(lrcBegin(), lrcEnd(), rhs.lrcBegin(), rhs.lrcEnd());
    }
    bool operator>(const BigUInt& rhs) const { return rhs < *this; }
    bool operator<=(const BigUInt& rhs) const;
    bool operator>=(const BigUInt& rhs) const { return !(*this < rhs); }

    /** Friends **/
    friend BigUInt operator+(size_t lhs, const BigUInt& rhs) { return rhs + lhs; }
    friend BigUInt operator-(size_t lhs, const BigUInt& rhs) { return BigUInt(lhs) - rhs; }
    friend BigUInt operator*(size_t lhs, const BigUInt& rhs) { return rhs * lhs; }
    friend BigUInt power(const BigUInt& base, size_t exponent);

    /***************** Builders *****************/
    static BigUInt        createRandom(size_t numberOfDigits);
    inline static BigUInt createRandomFromDecimalDigits(size_t orderOfMagnitude) {
        return createRandom(static_cast<size_t>(orderOfMagnitude * (std::log(10) / log(s_base)) + 1ul));
    }
    static inline BigUInt createWithRoom(size_t digitCount) {
        BigUInt result;
        result.reserve(digitCount);
        return result;
    }

    /***************** Output *****************/
    [[nodiscard]] std::string toString() const;
    [[nodiscard]] std::string toDecimalString() const;
    [[nodiscard]] std::string toBinaryString() const;
    friend std::ostream&      operator<<(std::ostream& os, const BigUInt& anInt);

private:
    /***************** Internal *****************/
    inline void init(size_t val) {
        if (val & s_highBits) {
            m_digits = std::vector<size_t>{val & s_lowBits, divideByBase(val)};
        } else {
            m_digits = std::vector<size_t>{val};
        }
    }
    void multiplyBySingleDigit(const size_t digit);
    void bubble(size_t startIndex = 0ul);
    void divideByLessThanBase(size_t factor);
    void square();

    [[nodiscard]] inline bool isZero() const {
        // assert(isWellFormed());
        return mostSignificantDigit() == 0ul;
    }
    [[nodiscard]] inline size_t value() const {
        return (digitCount() > 1ul) ? (leastSignificantDigit() + (mostSignificantDigit() << s_bitsPerDigit)) : leastSignificantDigit();
    }
    inline void reduceSizeByOneIfNeeded() {
        if (mostSignificantDigit() == 0ul) { resize(digitCount() - 1ul); }
        // assert(digitCount() >= 1ul);
    }
    inline void reduceSizeByTwoIfNeeded() {
        if (mostSignificantDigit() == 0ul) {
            resize(digitCount() - 1ul);
            reduceSizeByOneIfNeeded();
        }
    }

    /***************** Static helpers *****************/
    /** Vector **/
    static inline std::tuple<BigUInt, BigUInt, BigUInt, BigUInt> splitFour(rlcIterator begin, rlcIterator end, size_t i) {
        return {BigUInt({begin, begin + i}, false),
                BigUInt({begin + i, begin + 2ul * i}, false),
                BigUInt({begin + 2ul * i, begin + 3ul * i}, false),
                BigUInt({begin + 3ul * i, end}, false)};
    }
    static inline std::tuple<BigUInt, BigUInt, BigUInt> splitThree(rlcIterator begin, rlcIterator end, size_t i) {
        return {BigUInt({begin, begin + i}, false), BigUInt({begin + i, begin + 2ul * i}, false), BigUInt({begin + 2ul * i, end}, false)};
    }

    /** Addition **/
    static void carryAdditionViaIterators(rlIterator resultIt, size_t carry);
    static void carryUnitAdditionViaIterators(rlIterator resultIt);
    static void addViaIterators(rlIterator resultIt, rlcIterator rhsIt, rlcIterator rhsEnd);
    static void addMultipleViaIterators(rlIterator resultIt, rlcIterator rhsIt, rlcIterator rhsEnd, size_t multiplier);
    static void subtractViaIterators(rlIterator thisIt, rlcIterator rhsIt, rlcIterator rhsEnd);

    /** Multiplication **/
    static BigUInt multiply(const BigUInt& smaller, const BigUInt& larger);
    static void    multiplyBySingleDigitViaIterators(rlIterator resultIt, const rlIterator resultEnd, const size_t rhs);
    static void    karatsubaMultiplyViaIterators(
           rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd);
    static void splitOneMultiplicationViaIterators(
        rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd);
    inline static void
    multiplyViaIterators(rlIterator resultIt, rlcIterator lhsIt, const rlcIterator lhsEnd, rlcIterator rhsIt, const rlcIterator rhsEnd) {
        if ((rhsEnd - rhsIt) > lhsEnd - lhsIt) {
            multiplySortedViaIterators(resultIt, lhsIt, lhsEnd, rhsIt, rhsEnd);
        } else {
            multiplySortedViaIterators(resultIt, rhsIt, rhsEnd, lhsIt, lhsEnd);
        }
    }
    static void multiplySortedViaIterators(
        rlIterator resultIt, rlcIterator smallIt, const rlcIterator smallEnd, rlcIterator largeIt, const rlcIterator largeEnd);
    static void toomCook_3(rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd);
    static void toomCook_4(rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, rlcIterator largeEnd);
    static void schoolMultiply(rlIterator resultIt, rlcIterator smallIt, rlcIterator smallEnd, rlcIterator largeIt, size_t largeSize);

    /** Division **/
    static size_t  divisionSubRoutine(const lrcIterator leftToRightConstIt,
                                      const lrcIterator leftToRightConstEnd,
                                      const rlIterator  rightToLeftIt,
                                      const BigUInt&    divisor);
    static BigUInt longDivisionAfterAdjustingDivisor(BigUInt& dividend, const BigUInt& divisor);
    /** Comparison **/
    static bool lessThanShiftedRhsViaIterators(
        lrcIterator thisIt, const lrcIterator thisEnd, lrcIterator rhsIt, const lrcIterator rhsEnd, size_t trailingZeroesOfRhs);
    static bool
    lessThanViaIterators(const lrcIterator thisIt, const lrcIterator thisEnd, const lrcIterator rhsIt, const lrcIterator rhsEnd) {
        return lessThanShiftedRhsViaIterators(thisIt, thisEnd, rhsIt, rhsEnd, 0ul);
    }
    static bool
    greaterThanViaIterators(const lrcIterator thisIt, const lrcIterator thisEnd, const lrcIterator rhsIt, const lrcIterator rhsEnd) {
        return lessThanViaIterators(rhsIt, rhsEnd, thisIt, thisEnd);
    }
};
} // namespace big

#endif // __BIG_U_INT__H__