#include "timer.H"

MyTimer::MyTimer( )
{
    running_total = std::chrono::duration<double>::zero();
}

void MyTimer::start( )
{
    start_time = std::chrono::high_resolution_clock::now();
}

void MyTimer::end( )
{
    auto end_time = std::chrono::high_resolution_clock::now();
    duration = end_time - start_time;
}

void MyTimer::add( )
{
    running_total += duration;
}

double  MyTimer::get_total( )
{
    return running_total.count();
}