
#include "structure.H"

void start_clock( struct domain * theDomain ){
   theDomain->Wallt_init = time(NULL);
}
 
int count_cells( struct domain * theDomain ){

   int Nr = theDomain->Nr;

   int imin = 0;
   int imax = Nr;

   int Nc = imax-imin;

   return(Nc);
}

void generate_log( struct domain * theDomain ){
   time_t endtime = time(NULL);
   int seconds = (int) (endtime - theDomain->Wallt_init);
   
   int Nc = count_cells( theDomain );
   int Nt = theDomain->count_steps;

   double avgdt = (double)seconds/2./(double)Nc/(double)Nt;

   char logfile_filename[256] = "";
   strcat(logfile_filename, theDomain->output_prefix);
   strcat(logfile_filename, "times.log");

   FILE * logfile = fopen(logfile_filename,"w");
   fprintf(logfile,"Total time = %d sec\n",seconds);
   fprintf(logfile,"Number of cells = %d\n",Nc);
   fprintf(logfile,"Number of timesteps = %d (x%d)\n",Nt,2);
   fprintf(logfile,"Megazones per second = %.2e\n",1./(avgdt*1e6));
   fprintf(logfile,"Megazones per CPU second = %.2e\n",1./(avgdt*1e6));
   fprintf(logfile,"Time/zone/step = %.2e microseconds\n",(avgdt*1e6));
   fclose(logfile);

}

/*
void profiler_start( clock_t * prevtime , clock_t * currtime ){
   *prevtime = clock();
   *currtime = clock();
   printf("\n***\nProfiler Running...\n");
}

void profiler_report( char * event , clock_t * prevtime , clock_t * currtime ){
   *currtime = clock();
   printf("%s:\t\t%d ticks\n",event,(int)(*currtime-*prevtime) );
   *prevtime = *currtime;
}

void profiler_end(void){
   printf("***\n\n");
}
*/
