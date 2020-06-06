#ifndef __BIG_INT__H__
#define __BIG_INT__H__

#include "BigUInt.h"

#include <ostream>

namespace big {

    class BigInt {

    public:
        /***************** Constructors *****************/
        BigInt();
        BigInt(const size_t val);
        BigInt(const long long val);
        BigInt(const BigInt &other);
        explicit BigInt(const std::string &val);
        explicit BigInt(const BigUInt &other);
        BigInt(BigInt &&other) noexcept;
        explicit BigInt(BigUInt &&magnitude);

        /***************** Operators *****************/
        BigInt &operator=(size_t rhs);
        BigInt &operator=(long long rhs);
        BigInt &operator=(const BigUInt &rhs);
        BigInt &operator=(const BigInt &rhs);
        BigInt &operator=(BigInt &&rhs);
        BigInt &operator=(BigUInt &&rhs);

        /** Addition **/
        BigInt &operator+=(size_t rhs);
        BigInt &operator+=(long long rhs);
        BigInt &operator+=(const BigInt &rhs);
        BigInt &operator+=(const BigUInt &rhs);
        BigInt  operator+(size_t rhs) const;
        BigInt  operator+(const BigInt &rhs) const;

        /** Subtraction **/
        BigInt &operator-=(const BigInt &rhs);
        BigInt &operator-=(const BigUInt &rhs);
        BigInt  operator-(const BigInt &rhs) const;
        BigInt  operator-(const BigUInt &rhs) const;

        /** Multiplication **/
        BigInt &operator*=(size_t rhs);
        BigInt &operator*=(const BigInt &rhs);
        BigInt  operator*(size_t rhs) const;
        BigInt  operator*(const BigInt &rhs) const;

        /** Division **/
        BigInt &operator/=(size_t divisor);
        BigInt &operator/=(const BigInt &divisor);
        BigInt &operator/=(const BigUInt &divisor);
        BigInt  operator/(const BigInt &divisor) const;
        BigInt  operator/(const BigUInt &divisor) const;
        BigInt  operator/(size_t divisor) const;

        /** Comparison **/
        bool operator==(const BigInt &rhs) const { return m_magnitude == rhs.m_magnitude && m_isNegative == rhs.m_isNegative; }
        bool operator!=(const BigInt &rhs) const { return !(rhs == *this); }
        bool operator<(const BigInt &rhs) const;
        bool operator>(const BigInt &rhs) const { return rhs < *this; }
        bool operator<=(const BigInt &rhs) const;
        bool operator>=(const BigInt &rhs) const { return !(*this < rhs); }

        /** Friends **/
        friend BigInt operator+(size_t lhs, const BigInt &rhs);
        friend BigInt operator-(size_t lhs, const BigInt &rhs);
        friend BigInt operator*(size_t lhs, const BigInt &rhs);
        friend BigInt operator*(long long lhs, const BigInt &rhs);

        /***************** Output *****************/
        std::string          toString() const { return m_isNegative ? "=" : "" + m_magnitude.toString(); }
        friend std::ostream &operator<<(std::ostream &os, const BigInt &anInt);

        /***************** Vector functions *****************/
        void resize(size_t size);
        void reserve(size_t size);

        /***************** Accessors *****************/
        friend BigUInt;
        const BigUInt &magnitude() const { return m_magnitude; }

    private:
        void     negate() { m_isNegative = !m_isNegative; }
        BigUInt &magnitude() { return m_magnitude; }

        /***************** Data members *****************/
        BigUInt m_magnitude;
        bool    m_isNegative = false;
    };
} // namespace big
#endif // __BIG_INT__H__
