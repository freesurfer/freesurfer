#pragma once

#include <random>
#include <iostream>
#include <chrono>

class RandomSupply {
public:
  RandomSupply( unsigned int seed ) : generator(seed), mySeed(seed) {
  }

  RandomSupply() : generator() {
    this->mySeed = std::chrono::system_clock::now().time_since_epoch().count();
    this->generator.seed(this->mySeed);
  }

  template<typename T>
  T GetInteger(const T maxValue ) {
    auto randVal = this->generator();

    return (randVal % maxValue);
  }

  unsigned int GetSeed() const {
    return this->mySeed;
  }

  unsigned int Reseed() {
    return this->Reseed(std::chrono::system_clock::now().time_since_epoch().count());
  }

  unsigned int Reseed( unsigned int seed ) {
    this->generator.seed(seed);
    this->mySeed = seed;

    return this->mySeed;
  }
  
private:
  std::minstd_rand0 generator;
  unsigned int mySeed;
};
