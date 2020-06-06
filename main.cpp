#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

int main() {
    BenchMark::run({{2, 100},{4, 100},{10, 100}, {100, 100}, {1000, 100}, {10000, 50}, {100000, 10}, {1000000, 5}});
}
