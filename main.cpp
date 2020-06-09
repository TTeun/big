#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

using namespace big;

bool isPrime(const BigUInt& b) {
    for (BigUInt c = 3; c * c <= b; c += 2) {
        if (b % c == 0) { return false; }
    }
    return true;
}
int main() {
    BenchMark::run({{4, 400}, {10, 400}, {100, 400}, {1000, 400}, {10000, 50}});
}
