#pragma once

#include <chrono>
#include <stdexcept>

namespace kvl {
  class Stopwatch {
  public:
    Stopwatch() : running(false),
		  startCount(0),
		  start(),
		  elapsed() {}

    typedef long long CountType;
    typedef std::chrono::duration<CountType,std::micro> TimeSpan;

    void Start() {
      if( !running ) {
	this->start = std::chrono::high_resolution_clock::now();
	this->startCount++;
	this->running = true;
      } else {
	throw std::runtime_error("Attempted to start running Stopwatch");
      }
    }

    void Stop() {
      if( running ) {
	auto t2 = std::chrono::high_resolution_clock::now();
	this->elapsed += std::chrono::duration_cast<TimeSpan>(t2 - this->start);
	this->running = false;
      } else {
	throw std::runtime_error("Attempted to stop halted Stopwatch");
      }
    }

    template<typename Period = std::milli,typename Representation = float>
    Representation GetElapsedTime() const {
      typedef std::chrono::duration<Representation,Period> Target;
      auto t = std::chrono::duration_cast<Target>(this->elapsed);
      return t.count();
    }

    
    template<typename Period = std::milli, typename Representation = float>
    Representation GetAverageElapsedTime() const {
      typedef std::chrono::duration<Representation,Period> Target;
      auto t = std::chrono::duration_cast<Target>(this->elapsed);
      return t.count() / startCount;
    }

    void Reset() {
      this->running = false;
      this->startCount = 0;
      this->elapsed = TimeSpan::zero();
    }

  private:
    bool running;
    unsigned int startCount;
    std::chrono::high_resolution_clock::time_point start;
    TimeSpan elapsed;
  };
}
