#include "timer.h"


/// Returns the elapsed time in nanoseconds.
long Timer::nanoseconds() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() - begin).count();
}

/// Returns the elapsed time in milliseconds.
long Timer::milliseconds() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(clock::now() - begin).count();
}

/// Returns the elapsed time in seconds.
double Timer::seconds() {
    return std::chrono::duration_cast<std::chrono::duration<double>>(clock::now() - begin).count();
}

/// Returns the elapsed time in minutes.
double Timer::minutes() {
    return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<60>>>(clock::now() - begin).count();
}

/// Returns the elapsed time in minutes.
double Timer::hours() {
    return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<3600>>>(clock::now() - begin).count();
}

/// Resets the clock.
void Timer::reset() {
    begin = clock::now();
}


/// Returns the current date and time in the standard format.
///
/// For debugging and comparison purposes, it might be appropriate to always print the
/// same datetime string. Setting the 'FREESURFER_REPLACEMENT_FOR_CREATION_TIME_STRING'
/// environment variable to the desired replacement will override the true datetime
/// unless the `allowOverride` argument is set to `false`.

std::string currentDateTime(bool allowOverride) {
    const char* overrideString = getenv("FREESURFER_REPLACEMENT_FOR_CREATION_TIME_STRING");
    if (allowOverride && overrideString) {
        return overrideString;
    } else {
        char datetime[1000];
        time_t t = std::time(nullptr);
        struct tm *now = localtime(&t);
        std::strftime(datetime, 1000, "%c", now);
        return datetime;
    }
}
