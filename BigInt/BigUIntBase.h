#ifndef __DIGIT_VECTOR__H__
#define __DIGIT_VECTOR__H__

#include <cassert>
#include <ostream>
#include <vector>

typedef std::vector<size_t>::reverse_iterator       lrIterator;
typedef std::vector<size_t>::const_reverse_iterator lrcIterator;
typedef std::vector<size_t>::iterator               rlIterator;
typedef std::vector<size_t>::const_iterator         rlcIterator;

namespace big {

class BigUIntBase {
public:
    static const size_t s_bitsPerDigit    = std::numeric_limits<size_t>::digits / 2ul;
    static const size_t s_base            = 1ul << s_bitsPerDigit;
    static const size_t s_maxDigit        = s_base - 1ul;
    static const size_t s_lowBits         = (1ul << s_bitsPerDigit) - 1ul;
    static const size_t s_highBits        = std::numeric_limits<size_t>::max() ^ s_lowBits;
    static const size_t s_additionRoom    = std::numeric_limits<size_t>::max() - s_base;
    static const size_t s_decimalsInDigit = 9ul;
    static const size_t s_tenPowerInDigit = 1000000000ul;

    friend class BigInt;

    /***************** Utility members *****************/
    [[nodiscard]] static inline size_t divideByBase(size_t val) { return (val & s_highBits) >> s_bitsPerDigit; }

protected:
    /***************** Constructors *****************/
    BigUIntBase() = default;
    explicit BigUIntBase(std::vector<size_t>&& digits) : m_digits(std::move(digits)) { }

    /***************** Const functions *****************/
    [[nodiscard]] inline size_t digitAt(size_t index) const {
        // assert(index < m_digits.size());
        return m_digits.at(index);
    }
    [[nodiscard]] inline size_t leastSignificantDigit() const { return m_digits.front(); }
    [[nodiscard]] bool          isWellFormed() const;
    [[nodiscard]] inline size_t digitCount() const { return m_digits.size(); }
    [[nodiscard]] inline size_t mostSignificantDigit() const { return m_digits.back(); }
    [[nodiscard]] inline bool   isCorrectlySized() const {
        // assert(digitCount() > 0);
        return mostSignificantDigit() == 0ul ? (digitCount() == 1ul) : true;
    }

    /***************** Iterators *****************/
    [[nodiscard]] inline lrIterator  lrBegin() { return m_digits.rbegin(); }
    [[nodiscard]] inline lrIterator  lrEnd() { return m_digits.rend(); }
    [[nodiscard]] inline rlIterator  rlBegin() { return m_digits.begin(); }
    [[nodiscard]] inline rlIterator  rlEnd() { return m_digits.end(); }
    [[nodiscard]] inline lrcIterator lrcBegin() const { return m_digits.crbegin(); }
    [[nodiscard]] inline lrcIterator lrcEnd() const { return m_digits.crend(); }
    [[nodiscard]] inline rlcIterator rlcBegin() const { return m_digits.cbegin(); }
    [[nodiscard]] inline rlcIterator rlcEnd() const { return m_digits.cend(); }

    /***************** Vector functions *****************/
    inline void reserve(size_t size) { m_digits.reserve(size); }
    inline void resizeToFit() { resizeToFitVector(m_digits); }
    inline void resize(size_t size) { m_digits.resize(size); }
    static void resizeToFitVector(std::vector<size_t>& digits);
    static void copyViaIterators(rlcIterator begin, rlcIterator end, rlIterator resultIt);

    /***************** Data members *****************/
    std::vector<size_t> m_digits;
};
} // namespace big

#endif // __DIGIT_VECTOR__H__