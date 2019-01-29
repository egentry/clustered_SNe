#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>

class MyTimer
{
    public:
        MyTimer( );
        void start(  );
        void end( );
        void add( );
        double get_total( );

    protected:

    private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

    std::chrono::duration<double> duration; 
    std::chrono::duration<double> running_total;

};

#endif
