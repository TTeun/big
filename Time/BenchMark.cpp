#include "BenchMark.h"

#include "../BigInt/BigUInt.h"
#include "timer.h"

#include <gmpxx.h>
#include <iomanip>
#include <iostream>

using namespace big;

double BenchMark::multiply(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toBinaryString().c_str(), 2);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toBinaryString().c_str(), 2);
        {
            Timer t;
            mpz_mul(n1, n1, n2);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b *= a;
            bigTime += t.elapsed();
        }
        mpz_t n3;
        mpz_init(n3);
        mpz_set_str(n3, b.toBinaryString().c_str(), 2);
        assert(mpz_cmp(n1, n3) == 0);
    }
    return bigTime / gmpTime;
}

double BenchMark::add(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toBinaryString().c_str(), 2);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toBinaryString().c_str(), 2);
        {
            Timer t;
            mpz_add(n1, n1, n2);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b += a;
            bigTime += t.elapsed();
        }
        mpz_t n3;
        mpz_init(n3);
        mpz_set_str(n3, b.toBinaryString().c_str(), 2);
        assert(mpz_cmp(n1, n3) == 0);
    }
    return bigTime / gmpTime;
}

double BenchMark::divide(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen / 2ul);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toBinaryString().c_str(), 2);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toBinaryString().c_str(), 2);
        {
            Timer t;
            mpz_div(n1, n2, n1);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b /= a;
            bigTime += t.elapsed();
        }
        mpz_t n3;
        mpz_init(n3);
        mpz_set_str(n3, b.toBinaryString().c_str(), 2);
        assert(mpz_cmp(n1, n3) == 0);
    }
    return bigTime / gmpTime;
}

double BenchMark::modulo(size_t powerOfTen, size_t numberOfRepetitions) {
    double gmpTime = 0.0;
    double bigTime = 0.0;
    for (size_t i = 0; i != numberOfRepetitions; ++i) {
        const BigUInt a = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        BigUInt       b = BigUInt::createRandomFromDecimalDigits(powerOfTen);
        mpz_t         n1;
        mpz_init(n1);
        mpz_set_str(n1, a.toBinaryString().c_str(), 2);
        mpz_t n2;
        mpz_init(n2);
        mpz_set_str(n2, b.toBinaryString().c_str(), 2);
        {
            Timer t;
            mpz_mod(n1, n2, n1);
            gmpTime += t.elapsed();
        }
        {
            Timer t;
            b %= a;
            bigTime += t.elapsed();
        }
        mpz_t n3;
        mpz_init(n3);
        mpz_set_str(n3, b.toBinaryString().c_str(), 2);
        assert(mpz_cmp(n1, n3) == 0);
    }
    return bigTime / gmpTime;
}

void BenchMark::run(const std::vector<std::pair<size_t, size_t>> &digitCountAndRepetitions) {
    std::cout << "Decimal digits:\t";
    for (auto it : digitCountAndRepetitions) {
        std::cout << std::setfill(' ') << std::setw(15) << it.first;
    }
    std::cout << '\n';

    std::cout << "Add:\t\t\t";
    for (auto it : digitCountAndRepetitions) {
        std::cout << std::setfill(' ') << std::setw(15) << add(it.first, it.second);
    }
    std::cout << '\n';

    std::cout << "Multiply:\t\t";
    for (auto it : digitCountAndRepetitions) {
        std::cout << std::setfill(' ') << std::setw(15) << multiply(it.first, it.second);
    }
    std::cout << '\n';

    std::cout << "Divide:\t\t\t";
    for (auto it : digitCountAndRepetitions) {
        std::cout << std::setfill(' ') << std::setw(15) << divide(it.first, it.second);
    }
    std::cout << '\n';

    std::cout << "Modulo:\t\t\t";
    for (auto it : digitCountAndRepetitions) {
        std::cout << std::setfill(' ') << std::setw(15) << modulo(it.first, it.second);
    }
    std::cout << '\n';
}