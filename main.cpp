#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

int main() {
    BenchMark::run({{2,    400},
                    {4,    400},
                    {10,   400},
                    {100,  400},
                    {1000, 400}});
}
