#ifndef __DIGIT_VECTOR__H__
#define __DIGIT_VECTOR__H__

#include <cassert>
#include <ostream>
#include <vector>

typedef std::vector<size_t>::reverse_iterator lrIterator;
typedef std::vector<size_t>::const_reverse_iterator lrcIterator;
typedef std::vector<size_t>::iterator rlIterator;
typedef std::vector<size_t>::const_iterator rlcIterator;

namespace big {

    class BigUIntBase {

    public:
        static const size_t s_bitsPerDigit = std::numeric_limits<size_t>::digits / 2ul;
        static const size_t s_base = 1ul << s_bitsPerDigit;
        static const size_t s_maxDigit = s_base - 1ul;
        static const size_t s_additionRoom = std::numeric_limits<size_t>::max() - s_base;

        friend class BigInt;

    public:
        /***************** Constructors *****************/
        BigUIntBase() = default;

        explicit BigUIntBase(std::vector<size_t> &&digits) : m_digits(std::move(digits)) {}

        /***************** Const functions *****************/
        size_t digitAt(size_t index) const {
            assert(index < m_digits.size());
            return m_digits.at(index);
        }

        size_t leastSignificantDigit() const { return m_digits.front(); }

        bool isWellFormed() const;

        size_t digitCount() const { return m_digits.size(); }

        size_t mostSignificantDigit() const { return m_digits.back(); }

        bool isCorrectlySized() const;

        /***************** Iterators *****************/
        lrIterator lrBegin() { return m_digits.rbegin(); }

        lrIterator lrEnd() { return m_digits.rend(); }

        rlIterator rlBegin() { return m_digits.begin(); }

        rlIterator rlEnd() { return m_digits.end(); }

        lrcIterator lrcBegin() const { return m_digits.crbegin(); }

        lrcIterator lrcEnd() const { return m_digits.crend(); }

        rlcIterator rlcBegin() const { return m_digits.cbegin(); }

        rlcIterator rlcEnd() const { return m_digits.cend(); }

        /***************** Vector functions *****************/
        void reserve(size_t size) { m_digits.reserve(size); }

        void resizeToFit() { resizeToFitVector(m_digits); }

        void resize(size_t size) { m_digits.resize(size); }

        static void resizeToFitVector(std::vector<size_t> &digits);

        std::vector<size_t> shiftedCopy(size_t shiftAmount) const;

        static void copy(rlcIterator begin, rlcIterator end, rlIterator resultIt) {
            for (; begin != end; ++begin, ++resultIt) {
                *resultIt = *begin;
            }
        }

        /***************** Data members *****************/
        std::vector<size_t> m_digits;
    };
} // namespace big

#endif // __DIGIT_VECTOR__H__