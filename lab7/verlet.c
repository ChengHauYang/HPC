#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main() {
   // function declaration
   void Verlet(const int whichrun,
               const double u0,
               const double v0,
               const double dt,
               const int NumSteps);
   const double Tfinal = 12.0;
   const int NumSteps = 100;

   double dt=Tfinal/NumSteps;

   // run #1
   Verlet(1,1.0,0.4,dt,NumSteps);
   // run #2
   Verlet(2,1.0,0.6,dt,NumSteps);
   // run #3
   Verlet(3,1.0,0.8,dt,NumSteps);
   // run #4
   Verlet(4,-1.0,0.6,dt,NumSteps);
   // run #5
   Verlet(5,-1.0,0.4,dt,NumSteps);

return 0; }

void Verlet(const int whichrun,
             const double u0,
             const double v0,
             const double dt,
             const int NumSteps)
 {
    double u[NumSteps+1];
    double v[NumSteps+1];

    // Initial conditions
    u[0] = u0;
    v[0] = v0;

    double vstar;

    // Main loop
    for (int i=0;i<NumSteps;i++){
      vstar=v[i]+0.5*dt*(u[i]-pow(u[i],3));
      u[i+1]=u[i]+dt*vstar;
      v[i+1]=vstar+0.5*dt*(u[i+1]-pow(u[i+1],3));
    }

    // Print solution to file
    char filename[] = "outputX.data";
    filename[6] = whichrun+'0';

    // open file
    FILE* outfile = fopen(filename,"w");

    // output data
    for (int i=0;i<NumSteps+1;i++){
      fprintf(outfile,"%23.16e  ",u[i]);
      fprintf(outfile,"%23.16e\n",v[i]);
    }

    // close file
    fclose(outfile);

    // Call python script to plot
    system("python phase_plot.py");
}
