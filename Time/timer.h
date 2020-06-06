#ifndef __TIMER__H__
#define __TIMER__H__

#include <chrono>

class Timer {

public:
    Timer();
    ~Timer();
    double elapsed() const;

private:
    std::chrono::high_resolution_clock::time_point m_startTime = std::chrono::high_resolution_clock::now();
};

#endif // __TIMER__H__