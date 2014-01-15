#include <iostream>
#include <complex>

#include <bar.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>

#define CHARGE 1.0
// total no. of cycles
#define N_C 20.0
#define N_R 5.0
// absolute phase
#define PHI (M_PI/2.0)
// frequ. in au
#define W 0.17

// pulse duration
const double duration = N_C*2*M_PI/W;
// electric field amplitude, can be accurately changed
double E_0;

double vecpot_x(double time, int me);
double vecpot_y(double time, int me);
double vecpot_z(double time, int me);
double scalarpotx(double x, double y, double z, double time, int me);
double scalarpoty(double x, double y, double z, double time, int me);
double scalarpotz(double x, double y, double z, double time, int me);
double imagpot(long xindex, long yindex, long zindex, double time, grid g);
double field(double time, int me);

int main(int argc, char **argv)
{
  int lnointens, lintens, lno_of_ts, ts, me=0;
  double  intensity_min, intensity_max, dintensity, 
          intensity_cu, real_timestep, time, N, yield_N, yield_P,
          E_tot, gamma_K, I_p, U_p, z_expect;
  complex timestep;
  grid g, g_load;
  fluid ells;
  wavefunction staticpot, wf_load, wf, P;
  hamop hamilton;
  FILE *f1;

  lnointens = 11;
  intensity_min = 1e11;
  intensity_max = 1e12;

  g.set_dim(34);
  g.set_ngps(500, 5, 1); // 5 angular momenta 0,1,...,4
  g.set_delt(0.2);
  g.set_offs(0,0,0);

  g_load.set_dim(34);
  g_load.set_ngps(500, 1, 1); 
  g_load.set_delt(0.2,0.0,0.0);
  g_load.set_offs(0,0,0);

  dintensity = (intensity_max-intensity_min)/lnointens-1;

  for(lintens=0; lintens<lnointens; lintens++)
  {
    intensity_cu = intensity_min + dintensity*lintens;    

    // electric field strength, au for vecpot_*() in potentials.cc
    E_0 = sqrt(intensity_cu / 3.5095e16);

    ells.init(g.ngps_z());
    ells[0] = 0;

    hamilton.init(g,vecpot_x,vecpot_y,vecpot_z,
                    scalarpotx,scalarpoty,scalarpotz,
                    imagpot,field);

    staticpot.init(g.size()); 
    staticpot.calculate_staticpot(g,hamilton);

    wf.init(g.size()); 
    wf_load.init(g_load.size());

    f1 = fopen("hydrogen_1s.wf.txt","r");
    wf_load.init(g_load, f1, 0, 0);
    wf.regrid(g, g_load, wf_load);    
    fclose(f1);


    real_timestep = 0.2/4.0;
    lno_of_ts = int(duration/real_timestep) + 1;
    timestep=complex(real_timestep, 0.0);
    time = 0.0;

    // real-time-propagation loop
    for (ts=0; ts<lno_of_ts; ts++)
    {
      time = time + real_timestep;

      E_tot = real(wf.energy(0.0,g,hamilton,me,staticpot,CHARGE));
      P = wf.project(g, g_load, wf_load);
      N = wf.norm(g);
      z_expect = real(wf.expect_z(g));

      // if (ts%2000 == 0)
      //   printf("%15.10le %15.10le %15.10le %15.10le\n", time, E_tot,
      //     real(conj(P[0])*P[0]), N);

      wf.propagate(timestep, time, g, hamilton, me, staticpot, 0, CHARGE);
    }

    yield_N = (1.0 - N);
    yield_P = (1.0 - real(conj(P[0])*P[0]));

    printf("%15.10le %15.10le %15.10le\n", intensity_cu, yield_N, yield_P);
  }

  return 0;
}


double vecpot_x(double time, int me)
{
  return 0.0;
}  

double vecpot_y(double time, int me)
{
  return 0.0;
}

double vecpot_z(double time, int me)
{
  double result=0.0;
  double ww=W/(2.0*N_C);
 
  if ((time>0.0) && (time<duration))
  {
    result = E_0/(4.0*(W*W*W-4.0*W*ww*ww))*
    (-8.0*ww*ww*sin(PHI)- 2*(W*W-4*ww*ww)*sin(W*time+PHI)+W*((W+2*ww)*
    sin(W*time-2*ww*time+PHI)+(W-2*ww)*sin(W*time+2*ww*time+PHI)));
  }

  return result;
}

double scalarpotx(double x, double y, double z, double time, int me)
{
  return -CHARGE/x;
}

double scalarpoty(double x, double y, double z, double time, int me)
{ 
  return 0.0;
}

double scalarpotz(double x, double y, double z, double time, int me)
{
  return 0.0;
}

double imagpot(long xindex, long yindex, long zindex, double time, grid g)
{
  return 100*pow( (xindex + 0.5)/g.ngps_x(), 16);
} 

double field(double time, int me)
{
  return 0.0;
}

// g++-4.9 -O8 ionization.cc -o ground -I../../src -lqprop -lm -L../../build/lib


