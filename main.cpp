#include "BigInt/BigUInt.h"
#include "Time/BenchMark.h"

#include <iostream>

int main() {
    BenchMark::run({{2,       400},
                    {4,       400},
                    {10,      400},
                    {100,     400},
                    {1000,    400},
                    {10000,   50},
                    {100000,  10},
                    {1000000, 10}
                   });
}
//Add:			        2.87695         3.1717        2.33545        3.07758        5.03747        6.73031        6.97098        6.03473
