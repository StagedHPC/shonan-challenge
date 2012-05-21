#include <iostream>
#include <vector>

using namespace std;
const int N = 20;

/*
  Challenge : change the program02 to program3.
 */


/* 
   The initial simulator 
   written in parallel array computation notation
   that illustrates a typical finite-volume PDE solver.
*/
void program0(vector<double> input, vector<double>& output) {
  vector<double> input2(N);
  vector<double> wind(N);
  for (int i = 0; i < N; ++i) {
    input2[i]=input[i];
  }
  for (int i = 0; i < N-1; ++i) {
    wind[i] = input2[i] * input2[i+1];
  }
  for (int i = 1; i < N+1; ++i) {
    output[i] = input2[i] - wind[i] + wind[i-1];
  }
} /*   B/F = 16/3 
       This estimate assumes the bandwidth for intermediate storages
       are negligibly large. In reality, however, the the bandwidth
       for the intermediate storage of size N will possibly be close
       to that of the input and output, and estimate of B/F = 16/3
       is doubtful. 

       Also there's more simple problem the intermediates cost
       extra memory capacity.

       Can we remove the intermediate while keeping the B/F ?
  */

/* 
   The initial simulator, two timesteps fused into one
   written in parallel array computation notation
*/
void program02(vector<double> input, vector<double>& output) {
  vector<double> input2(N);
  vector<double> wind(N);
  vector<double> middleput(N);
  vector<double> wind_next(N);

  for (int i = 0; i < N; ++i) {
    input2[i]=input[i];
  }
  for (int i = 0; i < N-1; ++i) {
    wind[i] = input2[i] * input2[i+1];
  }
  for (int i = 1; i < N+1; ++i) {
    middleput[i] = input2[i] - wind[i] + wind[i-1];
  }
  for (int i = 0; i < N-1; ++i) {
    wind_next[i] = middleput[i] * middleput[i+1];
  }
  for (int i = 1; i < N+1; ++i) {
    output[i] = middleput[i] - wind_next[i] + wind_next[i-1];
  }
} /*   B/F = 16/6 
       You can theoretically improve B/F by factor of two by 
       using fusion, but again, this program requires intermediate 
       storages, and the improve in B/F is doubtful unless you
       can remove these intermediates.
  */


/*
  The initial simulator, one of the typical implementation
  that uses loop notation instead of parallel computation 
 */
void program1(vector<double> input, vector<double>& output) {
  for (int i = 1; i < N-1; ++i) {
    double a0 = input[i-1];
    double a1 = input[i];
    double a2 = input[i+1];
    double b0 = a0*a1;
    double b1 = a1*a2;
    double c0 = a1 - b1 + b0;
    output[i] = c0;
  }
} /*   B/F = 32/4 
       Worse than the theoretically best B/F because you couldn't
       sync using the intermediate storage
   */


/*
  The timestep-fused simulator in loop notation
 */
void program2(vector<double> input, vector<double>& output) {
  for (int i = 2; i < N-2; ++i) {
    double a0 = input[i-2];
    double a1 = input[i-1];
    double a2 = input[i];
    double a3 = input[i+1];
    double a4 = input[i+2];
    double b0 = a0*a1;
    double b1 = a1*a2;
    double b2 = a2*a3;
    double b3 = a3*a4;
    double c0 = a1 - b1 + b0;
    double c1 = a2 - b2 + b1;
    double c2 = a3 - b3 + b2;
    double d0 = c0 * c1;
    double d1 = c1 * c2;
    double e0 = c1 - d1 + d0;
    output[i] = e0;
  }
} /*   B/F = 48/14
       B/F has improved compared to program1,
       but compared to ideal program02, 
       this requires triple more bandwidth
       so it will possibly run three times slower.
  */

void program3(vector<double> input, vector<double>& output) {
  int i = 2;
  double a0 = input[i-2];
  double a1 = input[i-1];
  double a2 = input[i];
  double a3 = input[i+1];
  double a4 = input[i+2];
  double b0 = a0*a1;
  double b1 = a1*a2;
  double b2 = a2*a3;
  double b3 = a3*a4;
  double c0 = a1 - b1 + b0;
  double c1 = a2 - b2 + b1;
  double c2 = a3 - b3 + b2;
  double d0 = c0 * c1;
  double d1 = c1 * c2;
  double e0 = c1 - d1 + d0;
  output[i] = e0;
  for (int i = 3; i < N-1; ++i) {
    a3=a4; b2=b3; c1=c2; d0=d1;
    a4 = input[i+2];
    b3 = a3*a4;
    c2 = a3 - b3 + b2;
    d1 = c1 * c2;
    e0 = c1 - d1 + d0;
    output[i] = e0;
  }
} /*   B/F = 16/6 
       optimal, and also requires no intermediate data, so this 
       b/F estimate is quite realistic!
   */


int main () {
  // initial condition
  vector<double> initial(N);
  for(int i = 0; i < N; ++i) {
    initial[i] = 1.0 + 0.1 * i;
  }

  // double buffering
  vector<double> bufA(N), bufB(N);

  // does the result match?
  vector<double> final0, final02, final1, final2, final3;

  bufA = initial;
  program0(bufA, bufB);
  program0(bufB, bufA);
  final0 = bufA;

  bufA = initial;
  program02(bufA, bufB);
  final02 = bufB;


  bufA = initial;
  program1(bufA, bufB);
  program1(bufB, bufA);
  final1 = bufA;

  bufA = initial;
  program2(bufA, bufB);
  final2 = bufB;

  bufA = initial;
  program3(bufA, bufB);
  final3 = bufB;


  for(int i = 2; i < N-2; ++i) {
    cout << final0[i] 
         << "\t" << final02[i]
         << "\t" << final1[i]
         << "\t" << final2[i]
         << "\t" << final3[i]
         << endl;
  }
}


/*
  Although this particular challenge is too simple or not relevant
  to staging, I hope this example help guide you to general 
  context in HPC.

  HPC is starving of data. The main memory is too slow for 
  arithmetic units, and the cache is too small or missing.
  
  If we can capture the bulk of calculation
  that can be performed when the next data come,
  and automatically emit that portion as a program,
  it's great!!
*/

