#include <iostream>
#include <complex>

#include <bar.h>
#include <fluid.h>
#include <grid.h>
#include <hamop.h>
#include <wavefunction.h>
// #include <kbhit.h>
#include </Users/stas/code/qprop_ok/src/kbhit.cc>

double nuclear_charge = 1.0;  // hydrogen
double n_c=5.0;  // trapezoidal pulse: no. of cycles of constant intensity
                  // sinusoidal pulse: total no. of cycles
double n_r=0.0;   // trapezoidal pulse: up and down ramp over n_r cycles
double phi=M_PI / 2.0;  // absolute phase
double w=0.057;     // frequ. in au
double E_0=0.0;   // electric field amplitude 
                  // (will be assigned in hydrogen_re.cc (or ionization.cc))
double duration=n_c*2*M_PI/w;  // will be assigned when vecpot_z is called
double ampl_im=100.0;  // amplitude of imaginary absorbing potential

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
  // return vecpot_sin2_shape(time);  // or sin^2-shaped pulse

  double result=0.0;
  double ww;

  ww=w/(2.0*n_c);
  // duration=n_c*2*M_PI/w;
 
  if ((time>0.0) && (time<duration))
  {
    result=E_0/(4.0*(w*w*w-4.0*w*ww*ww))*(-8.0*ww*ww*sin(phi)-
    2*(w*w-4*ww*ww)*sin(w*time+phi)+w*((w+2*ww)*
    sin(w*time-2*ww*time+phi)+(w-2*ww)*sin(w*time+2*ww*time+phi)));
  }

  return result;
}

double scalarpotx(double x, double y, double z, double time, int me)
{  

  double result=-nuclear_charge/x; // -Z/r
  return result;
}

double scalarpoty(double x, double y, double z, double time, int me)
{
  double result=0.0;  
  return result;
}

double scalarpotz(double x, double y, double z, double time, int me)
{
  double result=0.0;
  return result;
}

double imagpot(long xindex, long yindex, long zindex, double time, grid g)
{
  if (ampl_im>1.0)
    return ampl_im*pow( (xindex + 0.5)/g.ngps_x(), 16);
  else
    return 0.0;
} 

double field(double time, int me)
{
  double result=0.0;
  return result;
}

int main(int argc, char **argv)
{
  double E_tot, curr_time, acc, E_tot_prev;
  complex timestep;
  int i, steps;
  grid g;
  fluid ells;
  hamop hamilton;
  wavefunction E_i, staticpot, wf;
  FILE *f1, *f2;

  g.set_dim(34); // 44 elliptical polariz., 34 linear polariz.
  g.set_ngps(500,1,1);  // <------------------------------------------- 
  g.set_delt(0.2/nuclear_charge);  // <-------------------------------- delta r
  g.set_offs(0,0,0);

  ells.init(g.ngps_z());
  ells[0] = 0; 

  E_i.init(g.ngps_z());

  // the Hamiltonian
  hamilton.init(g, vecpot_x, vecpot_y, vecpot_z, 
    scalarpotx, scalarpoty, scalarpotz,
    imagpot, field);

  // this is the linear and constant part of the Hamiltonian
  staticpot.init(g.size()); 
  staticpot.calculate_staticpot(g,hamilton);


  // *** new initialization ***
  // *** the wavefunction array 
  wf.init(g.size()); 
  wf.init(g, 1, 1.0, ells);
  wf.normalize(g);

  f1 = fopen("f1.txt","w");
  wf.dump_to_file_sh(g, f1, 1);
  fclose(f1);



  // ********************************************************
  // ***** imaginary propagation ****************************
  // ********************************************************

  timestep = complex(0.0, -1.0*g.delt_x()/4.0);  
  curr_time = 0.0;
  i = 0;
  steps = 640000;
  E_tot = 0.0;

  for (i; i<steps; i++)
  {
    curr_time += imag(timestep);
    E_tot_prev = E_tot;
    E_tot = real(wf.energy(0.0, g, hamilton, 0, staticpot, nuclear_charge));
    acc = fabs((E_tot_prev-E_tot)/(E_tot_prev+E_tot));
    wf.propagate(timestep, 0.0, g, hamilton, 0, staticpot, 0, nuclear_charge);
    wf.normalize(g);

    if (is_time(0.1))
      printf("% 20.15le % 20.15le %20.15le %ld\n", curr_time, E_tot, acc, i);

    if (kbhit(0))
      break;
  }

  f2 = fopen("f2.txt","w");
  wf.dump_to_file_sh(g, f2, 1);
  fclose(f2);

  return 0;
}

// g++-4.9 -O8 ground.cc -o ground -I../../src -lqprop -lm -L../../build/lib









