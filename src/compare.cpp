#include <iostream>
#include <math.h>
#include <R.h>

using namespace std;

extern "C" {

  void calcFDR(int* numFDR, double* stats, double* permStats, int* numGreater, int* numPermStats)
  {
    int counter = 0;
    int continueLoop = 1;
    int* isGreater =  NULL;
    isGreater = new int[numFDR[0]];

    for(int k = 0; k < numFDR[0]; k++)
      {
	isGreater[k] = 0;
      }

    for(int j = 0; j < numFDR[0]; j++)
      { 
	continueLoop = 1;
	while(continueLoop == 1)
	  {
	    continueLoop = 0;
	    if(counter < numPermStats[0])
	      {
		if(permStats[counter] > stats[j])
		  {
		    continueLoop = 1;
		    counter++;
		    isGreater[j] = isGreater[j] + 1;
		  }
	      }
	  }
      }

    numGreater[0] = isGreater[0];

    for(int i = 1; i < numFDR[0]; i++)
      { 
	numGreater[i] = isGreater[i] + numGreater[i-1];
      }
    delete [] isGreater;
  }


  void cumMax(float *stats, float *monotonicStats, int *length)
  {
    monotonicStats[length[0] - 1] = stats[length[0] - 1];
    int i;
    for(i=(length[0] - 2);i >= 0 ;i--){
      if(stats[i] > monotonicStats[i+1]){
      monotonicStats[i] = monotonicStats[i+1];
      }
    }
  }
}
