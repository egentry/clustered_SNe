
#ifndef _PROFILER_H_
#define _PROFILER_H_

#include "structure.H"

void start_clock(  struct domain * theDomain );
int  count_cells(  struct domain * theDomain );
void generate_log( struct domain * theDomain );


#endif