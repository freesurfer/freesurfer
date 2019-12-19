#include <boost/test/unit_test.hpp>

#include <thread>
#include <iostream>

#include "stopwatch.hpp"

BOOST_AUTO_TEST_SUITE( Stopwatch )

BOOST_AUTO_TEST_CASE( SimpleStartStop )
{
  kvl::Stopwatch timer;

  timer.Start();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.Stop();

  BOOST_CHECK_GE( timer.GetElapsedTime(), 1000 );
  BOOST_CHECK_EQUAL( timer.GetElapsedTime() * 1000, timer.GetElapsedTime<std::micro>() );
  BOOST_CHECK_EQUAL( timer.GetElapsedTime(), timer.GetAverageElapsedTime() );
  BOOST_CHECK_EQUAL( timer.GetAverageElapsedTime(), timer.GetAverageElapsedTime<std::micro>() / 1000 );
}

BOOST_AUTO_TEST_CASE( DoubleStartStop )
{
  kvl::Stopwatch timer;

  timer.Start();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.Stop();
  timer.Start();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.Stop();

  BOOST_CHECK_GE( timer.GetElapsedTime(), 2000 );
  BOOST_CHECK_GE( timer.GetAverageElapsedTime(), 1000 );
}

BOOST_AUTO_TEST_CASE( NoStartStarted )
{
  kvl::Stopwatch timer;

  timer.Start();
  BOOST_CHECK_THROW( timer.Start(), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( NoStopStopped )
{
  kvl::Stopwatch timer;

  BOOST_CHECK_THROW( timer.Stop(), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( Reset )
{
  kvl::Stopwatch timer;

  timer.Start();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.Stop();
  timer.Reset(); // Reset the timer
  timer.Start();
  std::this_thread::sleep_for(std::chrono::seconds(1));
  timer.Stop();

  BOOST_CHECK_GE( timer.GetElapsedTime(), 1000 );
  BOOST_CHECK_LE( timer.GetElapsedTime(), 1100 );
  BOOST_CHECK_GE( timer.GetAverageElapsedTime(), timer.GetElapsedTime() );
}

BOOST_AUTO_TEST_SUITE_END();
