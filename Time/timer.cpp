#include "timer.h"

Timer::Timer() {
}

Timer::~Timer() {
}

double Timer::elapsed() const {
    return std::chrono::duration_cast<std::chrono::duration<double>>(
               std::chrono::high_resolution_clock::now() - m_startTime)
        .count();
}
