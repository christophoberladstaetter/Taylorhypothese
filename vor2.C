// ===============================================================================
//                                      VOR2
// ===============================================================================
//
// 2D incompressible Navier-Stokes code in vorticity-streamfunction formulation
// with penalized Brinkman obstacle boundary conditions (for pipe flow)
// or periodic boundary conditions (for vortex merger; decaying turbulence; etc.).
//
// ( (c) Alexander Kendl, Innsbruck; vor2 version V1.0 of 01.04.2021)
//
// For instructions on compiling and running the code: see the 'vor2.readme' file
//
// ================================================================================

// include libraries
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fftw3.h>

// define global constants and variables:
const static double ZwPi = 2. * M_PI;
// nxmax, nymax: maximum nx and ny sizes, to initialise arrays:
const static int nxmax = 2 * 512 + 2, nymax = 512 + 2; // nx + 2 = 2^n + 2  ideal for fft
static bool bwall = true;
static int nx, ny, nxh, nyh, nx1, ny1, nx0, nin, ict, jct;
static int imbc, jmbc, ipbc, jpbc;
static int itstp, itmax, npar, iobs, isolv, incon, istart, isgn, iadvect, isorcnt = 0;
static double hy, hy2, hysq, hyd, xyz, dr, amp, wsrc, sigma, rrr, pen, vx, vy, ox, oy;
static double diff, dhyv, n00, lx, dt, dtt, ddtt;
static double aspect, eno, enwo, v0avg;
// fields for stream function phi, and vorticity www:
static double phi[nxmax][nymax], ph0[nxmax][nymax], www[nxmax][nymax];
static double w00[nxmax][nymax], obs[nxmax][nymax], wall[nxmax][nymax];
static double e[nxmax][nymax]; // Definiere Energie Matrix
// stored back values of vorticity:
static double ww0[nxmax][nymax], ww1[nxmax][nymax], ww2[nxmax][nymax];
static double dd0[nxmax][nymax], dd1[nxmax][nymax], dd2[nxmax][nymax];
// spatial derivative terms and back values in vorticity equation
static double fw0[nxmax][nymax], fw1[nxmax][nymax], fw2[nxmax][nymax];
// FFT kernels
static double cvor[nxmax][nymax];
// terms in SOR solver
static double uu0[nxmax][nymax];
// additional fields:
static double obs_core[nxmax][nymax]; // field for additional analysis of flow through the obstacla
// functions are first defined here, see further below main part:
void initpar(void);
void arakawa(double (*uuu)[nymax], double (*vvv)[nymax], double (*ww)[nymax]);
void upwind(double (*uuu)[nymax], double (*vvv)[nymax], double (*ww)[nymax]);
void simple(double (*uuu)[nymax], double (*vvv)[nymax], double (*ww)[nymax]);
void laplace(double (*ww)[nymax], double (*pp)[nymax], int ib, int io);
void sor2p(double (*qqq)[nymax], double (*uuu)[nymax]);
void mirror(double (*arrin)[nymax]);
void diagnose(double ttt, double t00, int it, int itmax,
              double (*pp)[nymax], double (*ww)[nymax], double (*dd)[nymax],int, double, double, double, double );

// --------------------------- Main programme part ---------------------------------
double enk_i;
double reynold;
int main(void)
{
  int i, j, l, m, n, ik, jk, it, is, itn, itend, i0, j0, i1, iz, jz, jz2,iout, jout;
  int im, ip, jm, jp, im1, ip1, im2, ip2, nk;
  double c0 = 18. / 11., c1 = 9. / 11., c2 = 2. / 11., cf = 6. / 11.;
  double ttt, t00, phase, zuf, zuf1, zuf2, pbar, visc, bob, dox, doy;
  double xi, yi, dist, width, turb;
  
  double arw[nxmax][nymax], hyvw[nxmax][nymax], visw[nxmax][nymax];
  double ard[nxmax][nymax], hyvd[nxmax][nymax], visd[nxmax][nymax];
  char s[80], str[80];
  char const *fn, *fw;
  FILE *f, *g, *h, *g1, *g2, *g3, *ft, *g4, *g5,*g7, *ene, *grid;
  int counter = 0; // this is my file counter

  // read input file and initialise some parameters:
  initpar(); // see subroutine specified below

  // initialise parallelised FFT:
  // (on FFTW3 C subroutine library: see 'https://www.fftw.org')
#ifdef _OPENMP
  fftw_plan_with_nthreads(npar);
#endif
  fftw_complex wwk[nx][ny], ppk[nx][ny], wwx[nx][ny], ppx[nx][ny];
  fftw_plan hinp, herp;
  hinp = fftw_plan_dft_2d(nx, ny, &wwx[0][0], &wwk[0][0], FFTW_FORWARD, FFTW_MEASURE);
  herp = fftw_plan_dft_2d(nx, ny, &ppk[0][0], &ppx[0][0], FFTW_BACKWARD, FFTW_MEASURE);

  // initialise more variables and parameters:

  printf("| initialisation ... ");
  t00 = 0.;

  if (incon == 0) // pipe flow b.c. parameters with in/out flow and walls
  { // incon... initial flow condition
    nx0 = nxh; //nxh = nx / 2 und nx... bei pipe-flow ist es das init grid in x*2 
    i0 = nxh;
    imbc = 0; // imboundarycond
    ipbc = nx1; // inoundary maximum nx1 = nx - 1
    jmbc = 0;
    jpbc = ny1; // ny1 = ny - 1
  }
  else // periodic b.c parameters in x and y
  {
    imbc = nx1;
    ipbc = 0;
    jmbc = ny1;
    jpbc = 0;
    nx0 = 0;
    i0 = 0;
  };

  iz = nxh + (nx1 - nxh) / 2;
  jz = nyh / 10.; // this defines the place of the grid
  jz2 = nyh / 2; // this defines the place of the grid
  i1 = nx1; // nx1 = nx - 1 wobei nx grid points in x und bei incon==0 --> nx = nx*2 

  // Define several flow initial conditions by stream function phi:
//-------------------------------------------------------
  if (incon == 0) // linear vorticity www(x) ~  parabolic pipe flow profile vy(x) :
  { 
    // hier wird von der Hälfte bis ganz nach oben in x iteriert
    for (i = nx0; i <= i1; ++i)
    {
      xi = i - iz - 0.5;
      for (j = 0; j <= ny1; ++j)
      {
        ph0[i][j] = amp * (-xi / double(nxh / 2) + xi * xi * xi / (3. * nxh * nxh * nxh / 8.));
      }
    }
  }
//---------------------------------------------------------
  if (incon == 1) // double vortex (vortex merger or co-advection)
  {
    width = sigma * (nx1 - nx0);
    dist = 2. * width;
    for (i = nx0; i <= nx1; i++)
      for (j = 0; j <= ny1; j++)
      {
        xi = i - (nx1 - nx0) / 2;
        yi = j - ny1 / 2 - dist;
        www[i][j] = amp * exp(-(xi * xi + yi * yi) / (width * width));
        yi = j - ny1 / 2 + dist;
        www[i][j] += isgn * amp * exp(-(xi * xi + yi * yi) / (width * width));
      }
  }

  if (incon == 2) // perturbation field for decaying turbulence
  {
    nk = 16;
    double phs[nk + 1][nk + 1];
    // srand(time(NULL)); // + getpid());
    for (ik = 1; ik <= nk; ik++)
      for (jk = 1; jk <= nk; jk++)
      {
        zuf = 0.02 * (rand() % 50);
        phs[ik][jk] = cos(ZwPi * zuf) * amp / sqrt(1. + pow((3.125 * (ik * ik + jk * jk) / (nk * nk)), 4));
      }
    zuf = 0.1 * (rand() % 10);
    zuf1 = 0.1 * (rand() % 10);
    zuf2 = 0.1 * (rand() % 10);
    for (i = nx0; i <= nx1; i++)
      for (j = 0; j <= ny1; j++)
      {
        for (ik = 1; ik <= nk; ik++)
          for (jk = 1; jk <= nk; jk++)
          {
            turb = phs[ik][jk] * cos(ZwPi * zuf) * sin(ZwPi * ik * (i / double(nx)) + ZwPi * zuf1) * sin(ZwPi * jk * (j / double(ny)) + ZwPi * zuf2);
            ph0[i][j] += amp * turb / double(nk);
          }
        // for testing:
        // ph0[i][j] = amp*sin(2.*ZwPi*(i-nx0)/double(nx1-nx0+1))*sin(2.*ZwPi*j/double(ny1+1));
      }
  }

  if (incon == 3) // Kelvin-Helmholtz instabilitystatic double obs_core[nxmax][nymax]; // field for additional analysis of flow through the obstacla
  {
    width = sigma * (nx1 - nx0) / 2;

    for (i = nx0; i <= nx1; i++)
      for (j = 0; j <= ny1; j++)
      {
        xi = i - (nx1 - nx0) / 2;
        ph0[i][j] = amp * exp(-xi * xi / (width * width));
        ph0[i][j] += 0.0000001 * sin(1. * ZwPi * j / double(ny1 + 1)) * exp(-xi * xi / (0.1 * width * width));
      }
  }

  if (incon == 4) // single vortex (to show spreading depending on Re)
  {
    width = sigma * (nx1 - nx0);
    for (i = nx0; i <= nx1; i++)
      for (j = 0; j <= ny1; j++)
      {
        xi = i - (nx1 - nx0) / 2;
        yi = j - ny1 / 2;
        www[i][j] = amp * exp(-(xi * xi + yi * yi) / (width * width));
      }
  }


//----------------- here nice stuff happens--------------


  // copy start value of phi to evolving field:
  for (i = nx0; i <= i1; ++i)
    for (j = 0; j <= ny1; ++j)
      phi[i][j] = ph0[i][j];

  if ((incon != 1) && (incon != 4))
  {
    // obstacle and wall b.c.
    // recall nx0 = nxh; //nxh = nx / 2 und 
    // i1 = nx1; // nx1 = nx - 1 wobei nx grid points in x und bei incon==0 --> nx = nx*2 
    // iz = nxh + (nx1 - nxh) / 2;
    for (i = nx0; i <= nx1; ++i) 
    {
      // xi = i - iz - 0.5;
      xi = i - iz - 1; // slight asymmetry helps to form vortex street
      for (j = 0; j <= ny1; ++j)
      {
        obs[i][j] = 0.;
        wall[i][j] = 0.;
        yi = j - jz;
        int yi2 = j-jz2;
        if (bwall)
        {
          if ((i < nxh + 2) || (i > nx1 - 2))
            wall[i][j] = 1.;
        }              // x pipe wall boundaries
        if (iobs == 1) // plate obstacle----------------------- Do some changes to obtain a nice "fence" which makes lots of turbulece 23.März Oberladstätter
        {
          if ((fabs(xi) < (sigma * nxh)) && (fabs(yi) <= 2.))
            obs[i][j] = 1.; // Changes have to be done here
          if (sigma == 1.)
          {
            obs[i][j] = 0.;
            if ( ((i-nx0) > 0.00*(nx1-nx0)) && ((i-nx0) < 0.10*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.12*(nx1-nx0)) && ((i-nx0) < 0.14*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.16*(nx1-nx0)) && ((i-nx0) < 0.18*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.20*(nx1-nx0)) && ((i-nx0) < 0.22*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.24*(nx1-nx0)) && ((i-nx0) < 0.26*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.28*(nx1-nx0)) && ((i-nx0) < 0.30*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.32*(nx1-nx0)) && ((i-nx0) < 0.34*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.36*(nx1-nx0)) && ((i-nx0) < 0.38*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.40*(nx1-nx0)) && ((i-nx0) < 0.42*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.44*(nx1-nx0)) && ((i-nx0) < 0.46*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.48*(nx1-nx0)) && ((i-nx0) < 0.50*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.52*(nx1-nx0)) && ((i-nx0) < 0.54*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.56*(nx1-nx0)) && ((i-nx0) < 0.58*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.60*(nx1-nx0)) && ((i-nx0) < 0.62*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.64*(nx1-nx0)) && ((i-nx0) < 0.66*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;
            if ( ((i-nx0) > 0.68*(nx1-nx0)) && ((i-nx0) < 0.70*(nx1-nx0)) && (fabs(yi)<=2.))
             obs[i][j] = 1.;
            if ( ((i-nx0) > 0.72*(nx1-nx0)) && ((i-nx0) < 0.74*(nx1-nx0)) && (fabs(yi)<=2.))
             obs[i][j] = 1.;
            if ( ((i-nx0) > 0.76*(nx1-nx0)) && ((i-nx0) < 0.78*(nx1-nx0)) && (fabs(yi)<=2.))
             obs[i][j] = 1.;
            if ( ((i-nx0) > 0.80*(nx1-nx0)) && ((i-nx0) < 0.82*(nx1-nx0)) && (fabs(yi)<=2.))
             obs[i][j] = 1.;
            if ( ((i-nx0) > 0.84*(nx1-nx0)) && ((i-nx0) < 0.86*(nx1-nx0)) && (fabs(yi)<=2.))
             obs[i][j] = 1.;
            if ( ((i-nx0) > 0.88*(nx1-nx0)) && ((i-nx0) < 0.90*(nx1-nx0)) && (fabs(yi)<=2.))
             obs[i][j] = 1.;
            
            if ( ((i-nx0) > 0.92*(nx1-nx0)) && ((i-nx0) < 1.0*(nx1-nx0)) && (fabs(yi)<=2.))
              obs[i][j] = 1.;

        
          }
        }
        if (iobs == 2)
        { // circular obstacle
          rrr = sqrt(xi * xi + yi * yi);
          if (rrr < (sigma * nxh))
            obs[i][j] = 1.;
        }
      }
    }

    // Print obstacle to file
    g7 = fopen("simdata1/obstacle.dat", "w");
    for (i=nx0; i<=nx1; ++i){
      for (j=0; j<=ny1; ++j){
        obs_core[i][j] = obs[i][j];
        // narrow edge (4 neiggbours define core)
        // if (obs[i][j] && obs[i+1][j] && obs[i-1][j] && 
        // obs[i][j+1] && obs[i][j-1]) obs_core[i][j] = 2.;
        // wider edge (8 neighbours define core)
        if (obs[i][j] && obs[i+1][j] && obs[i+1][j+1] && obs[i][j+1] && 
            obs[i-1][j+1] && obs[i-1][j] && obs[i-1][j-1] && obs[i][j-1] && 
            obs[i+1][j-1]) obs_core[i][j] = 2.;
        fprintf(g7, "%d  %d  %.1e\n", i, j, obs_core[i][j]);
      }
      fprintf(g7, "\n");
    }
    fclose(g7);


    // store initial ph0 = phi(t=0) and w00 = www(t=0) for later re-use
    if (incon == 0)
    {
      mirror(ph0); // "invisible x mirror domain" reduces asymmetry for good FFT in pipe flow
      mirror(phi);
    }
    // Man hat die phi--> flowstream parabloisch definiert
    // Nun Berechnung der vorticity aus der streamfunction
    // omega = laplace(phi)
    laplace(w00, phi, 0, 0); // get initial vorticity from initial flow stream function phi
    for (i = 0; i <= nx1; ++i)
      for (j = 0; j <= ny1; ++j)
      {
        w00[i][j] = -w00[i][j]; // minus by definition of streamfunction
        www[i][j] = w00[i][j]; // Here the w00 vorticity
      };
  }

  v0avg = 0.; // average of initial pipe flow profile: used to subtract in visualisation
  for (i = nx0 + 1; i <= nx1 - 1; i++)
    for (j = 0; j <= ny1; j++)
    {
      v0avg += -.5 * hy * (phi[i + 1][j] - phi[i - 1][j]) / double(ny * (nxh - 2));
    }

  // set artificial "previous" time values for first step (in multi-step solver)
  for (i = 0; i <= nx1; ++i)
    for (j = 0; j <= ny1; ++j)
    {
      ww0[i][j] = www[i][j];
      ww1[i][j] = www[i][j];
      ww2[i][j] = www[i][j];
      fw1[i][j] = 0.;
      fw2[i][j] = 0.;
    }

  eno = 0.; // energy diagnostics: helps to identify crashes or numerical instabilities
  for (i = 0; i <= nx1; i++)
    for (j = 0; j <= ny1; j++)
      eno += www[i][j] * www[i][j];
  eno *= .5 * xyz; //   xyz = 1. / double(nx * ny); --> look at initialize parameters

  // so you can kind a continue with the simulation from the end of the previous simulation
  if (istart == 2) // restart from saved end output of previous run
  {
    printf("| continued data set... \n");
    g = fopen("restart.dat", "r");
    fscanf(g, "%s ", str);
    t00 = atof(str);
    while (!feof(g))
    {
      fscanf(g, "%s ", str);
      i = atoi(str);
      fscanf(g, "%s ", str);
      j = atoi(str);
      fscanf(g, "%s ", str);
      eno = atof(str);
      fscanf(g, "%s ", str);
      phi[i][j] = atof(str);
      fscanf(g, "%s ", str);
      www[i][j] = atof(str);
      fscanf(g, "%s ", str);
      ww1[i][j] = atof(str);
      fscanf(g, "%s ", str);
      ww2[i][j] = atof(str);
      fscanf(g, "%s ", str);
      fw1[i][j] = atof(str);
      fscanf(g, "%s ", str);
      fw2[i][j] = atof(str);
    }
    fclose(g);
    for (i = 0; i <= nx1; i++)
      for (j = 0; j <= ny1; j++)
      {
        ww0[i][j] = www[i][j];
      }
  }

  // control output of initial profiles etc:
  diagnose(0., 0., it, itmax, phi, www, www,counter, reynold,dhyv, dt, itstp);

  // exit(1);

  // initialise time stepping:
  printf("| start time stepping ...\n");
  itn = 0;
  if (itstp != 0)
    itend = itmax / itstp;
  if (itstp == 0)
    itend == 0;
  dtt = dt * cf; // dt is defined in the inip file; cf = 6. / 11.;
  ddtt = double(itstp) * dt; // itspt inner steps between outputs

  if (itmax == 0)
    return (0);

  // Time steps  --------------------------------------------------------------
  // 3rd order Karniadakis time scheme;
  // Arakawa brackets or simple upwinding for advective term;
  // standard finite differencing for other terms

  // itend = itmax / itstp
  // itstp(steps between output)
  // dt (time step)
  // itmax (total steps)
  for (it = itn / itstp; it < itend; ++it) // total time loop (output step)
  { 
    for (is = 0; is < itstp; ++is) // inner time loop (without outputs)
    {
      ttt = t00 + (it + 0) * ddtt + (is + 1) * dt;

      // advection term: compute Poisson brackets [phi,vor]:
      // advection solver is specified with the iadvect condition
      // Hier lösen der Poissonklammer
      if (iadvect == 0)
        upwind(phi, www, arw); // simple low-order upwind for pipe flows
      if (iadvect == 1)
        simple(phi, www, arw); // simple low-order centered scheme
      if (iadvect == 2)
        arakawa(phi, www, arw); // energy and enstrophy conserving, centered

      // viscosity and hyperviscosity terms:
      if (diff != 0.)
        laplace(visw, www, 0, 0); // for Reynolds number term
      if (dhyv != 0.) // müsste hier nicht zweimal Laplace sein??
        laplace(hyvw, visw, 0, 0); // adds to numerical stability for turbulence

        // update the vorticity field www:
    #ifdef _OPENMPnx
    #pragma omp parallel for private(j, jm, jp, visc, bob, vx, vy, ox, oy, dox, doy)
    #endif
      for (i = 0; i <= nx1; i++) // iteriert über alle x
      {
        im = (i == 0) ? imbc : i - 1;
        ip = (i == nx1) ? ipbc : i + 1;
        for (j = 0; j <= ny1; j++) // iteriert über ganzes Gitter in y
        {
          jm = (j == 0) ? jmbc : j - 1;
          jp = (j == ny1) ? jpbc : j + 1;

          // "physical" ~ (1/Re) viscosity term
          visc = +diff * visw[i][j];

          // artificial hyperviscosity ~ k^4 for numerical stabilization
          visc += -dhyv * hyvw[i][j]; // (take dyhv as small as possible)

          // treat obstacles and walls by Brinkman penalization method
          // (walls: zero velocities but here with finite initial vorticity)
          // [See e.g.: Spietz et al., J. Comp. Phys. 336, 261 (1017), in eqs. (1)-(6)]
          bob = -pen * obs[i][j] * www[i][j];                 // obstacle penalty 1st term
          bob += -pen * wall[i][j] * (www[i][j] - w00[i][j]); // wall penalty with driving
          // bob+=  - pen*wall[i][j]*(www[i][j]); // wall penalty without drive
          if ((iobs != 0) && (obs[i][j] == 1.)) // obstacle penalty 2nd term
          {
            ox = 0.;
            dox = hy * (obs[ip][j] - obs[i][j]);
            if (dox < 0.)
            {
              ox = dox;
              vy = -hy * (phi[ip][j] - phi[i][j]);
            }
            dox = hy * (obs[i][j] - obs[im][j]);
            if (dox > 0.)
            {
              ox = dox;
              vy = -hy * (phi[i][j] - phi[im][j]);
            }
            oy = 0.;
            doy = hy * (obs[i][jp] - obs[i][j]);
            if (doy < 0.)
            {
              oy = doy;
              vx = +hy * (phi[i][jp] - phi[i][j]);
            }
            doy = hy * (obs[i][j] - obs[i][jm]);
            if (doy > 0.)
            {
              oy = doy;
              vx = +hy * (phi[i][j] - phi[i][jm]);
            }
            bob += -pen * (ox * vy - oy * vx);
          }
          if (wall[i][j] == 1.) // wall penalty: no 2nd term
          {
            // by definition at wall boundary: vx = 0.; vy = 0.;
            // and therefore bob += - pen*(ox*vy - oy*vx) == 0;
            bob += 0.;
          }

          // "right-hand side" of the time evolution equation:
          // arw ist die gelöste Poissonklammer
          fw0[i][j] = -arw[i][j] + visc + bob;

          // Time-step update: here we use a  higher-order splitting
          // "Karniadakis" Adams-Bashforth type multi-step method
          // which requires saving the www field of the last two time steps
          // [Reference: G.E. Karniadakis, J. Comp. Phys. 97, 414 (1991)]

          // Anschreiben der Navier stokes und eben das lösen
          // dtWWW= 1/Re laplace(WWW) - [Phi,WWW] 
          www[i][j] = c0 * ww0[i][j] - c1 * ww1[i][j] + c2 * ww2[i][j] + dtt * (3. * fw0[i][j] - 3. * fw1[i][j] + fw2[i][j]);
        }
      }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif

      if (incon == 0)
        mirror(www);

      // inflow and outflow b.c. for pipe flow
      if (incon == 0)
      {
        for (i = 0; i <= nx1; i++)
          for (j = 0; j <= 5; j++)
          {
            www[i][j] = w00[i][j]; // nimm die ersten 5 Felder von den inits--> Einfluss
          };
        for (i = 0; i <= nx1; i++)
          for (j = ny1 - 2; j <= ny1; j++)
          {
            www[i][j] = www[i][ny1 - 3]; // nimm die letzten zwei Felder und setzte sie mit dem drittvorletzen gleich
          };
      }

      // solve Poisson equation: ----------------------------

      if (isolv == 1)
      {
        // SOR iterative Poisson solver for phi from www (slow)
        // Omega = laplace(phi)
        sor2p(www, phi);
        pbar = 0.;
        for (j = 0; j <= ny1; ++j)
          pbar += phi[nx1][j];
        pbar /= (ny);
        for (i = 0; i <= nx1; ++i)
          for (j = 0; j <= ny1; ++j)
            phi[i][j] -= pbar;
      }
      else if (isolv == 0) // here the fouriertransform happens
      {
        // FFT Poisson solver in k-space for phi from www
        // (fast in parallelisation; needs fftw3 include)
        #ifdef _OPENMP
        #pragma omp parallel for private(j)
        #endif
        // just for easy reading --> nx1 = nx-1
        for (i = 0; i <= nx1; ++i) // über alle x iterieren
          for (j = 0; j <= ny1; ++j) // über alle y iterieren
            wwx[i][j][0] = www[i][j];
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        //  hinp = fftw_plan_dft_2d(nx, ny, &wwx[0][0], &wwk[0][0], FFTW_FORWARD, FFTW_MEASURE);
        //  herp = fftw_plan_dft_2d(nx, ny, &ppk[0][0], &ppx[0][0], FFTW_BACKWARD, FFTW_MEASURE);
        // Change from the wwx into the wwk space 
            fftw_execute(hinp); // do forward transform
        #ifdef _OPENMP
        #pragma omp parallel for private(j)
        #endif
        for (i = 0; i <= nx1; ++i)
          for (j = 0; j <= ny1; ++j)
          { // warum hier [0] und [1]
            // Berechne aus neuem vorticity www neues streamfunction phi
            ppk[i][j][0] = wwk[i][j][0] * cvor[i][j]; // entspricht dem k²?
            ppk[i][j][1] = wwk[i][j][1] * cvor[i][j];
          }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
        fftw_execute(herp); // do back transform
        #ifdef _OPENMP
        #pragma omp parallel for private(j)
        #endif
        for (i = 0; i <= nx1; ++i)
          for (j = 0; j <= ny1; ++j)
            phi[i][j] = ppx[i][j][0]; // neue streamfunction
        #ifdef _OPENMP
        #pragma omp barrier
        #endif
      }

      // minus by definition of stream function:
      for (i = 0; i <= nx1; i++)
        for (j = 0; j <= ny1; j++)
        {
          phi[i][j] *= -1.;
        };

          // multi-step memory: remember values two steps backwards
    #ifdef _OPENMP
    #pragma omp parallel for private(j)
    #endif
      for (i = 0; i <= nx1; ++i)
        for (j = 0; j <= ny1; ++j)
        {
          ww2[i][j] = ww1[i][j];
          ww1[i][j] = ww0[i][j];
          ww0[i][j] = www[i][j];
          fw2[i][j] = fw1[i][j];
          fw1[i][j] = fw0[i][j];
        }
    #ifdef _OPENMP
    #pragma omp barrier
    #endif

    } // ... end inner time loop ..............................................

    // low-time resolution diagnostics: energies and snapshot outputs
    // (output to file is slow, so only done in every n-th time step)
    // this is my counter for my outputfiles

    if(ttt>1)
      { counter=counter+1; 
        diagnose(ttt, t00, it, itmax, phi, www, www, counter, reynold, dhyv, dt, itstp);}

  } // ... end outer time loop ..........................................

  // End time step --------------------------------------------------------

  // write complete restart file
  g = fopen("restart.dat", "w");
  fprintf(g, "%5e  \n", ttt);
  for (i = 0; i <= nx1; i++)
    for (j = 0; j <= ny1; j++)
      fprintf(g, "%d  %d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  \n",
              i, j, eno, phi[i][j], www[i][j], ww1[i][j], ww2[i][j], fw1[i][j], fw2[i][j]);
  fclose(g);

#ifdef _OPENMP
  fftw_cleanup_threads();
#endif

  printf("| VOR2 end.\n\n");
}

// END OF MAIN ================================================================

// Specify functions:

// ---------------------------------------------------------------------------
// Read input data file "vor2.inp" and initialise some parameters
void initpar(void)
{
  int i, j, ith;
  char s[80];
  FILE *g, *f;
  double rel, para[50];

  for (int i = 0; i <= 40; i++)
    printf("_");
  printf("\n| VOR2 start...\n");

  // this sequentially reads all numbers appearing behind any "=" sign in vor2.inp
  // caution: be aware of changing the numbering in case of additions to vor2.inp
  g = fopen("vor2.inp", "r");
  i = 0;
  while (!feof(g))
  {
    fscanf(g, "%s ", s);
    if (strstr(s, "="))
      i++;
  };
  rewind(g);
  for (j = 1; j <= i; j++)
  {
    while (!strstr(s, "="))
      fscanf(g, "%s ", s);
    fscanf(g, "%s ", s);
    para[j] = atof(s);
  };
  fclose(g);

  if (para[1] == 0.)
    bwall = false;

  iobs = int(para[2]); // obstacle type

  istart = int(para[3]); // start condition: 0 = new, 2 = restart

  diff = 1. / para[4]; // inverse Reynolds number (viscosity term prefactor)

  reynold = para[4]; // read in the reynoldsnumber 

  if (para[4] == 0.)
    diff = 0.;
  if (iobs != 0)
    diff *= para[6]; // scale Re with obstacle size

  incon = int(fabs(para[5])); // initial flow condition:
                              // 0 : parabolic pipe flow profile
                              // 1 : dual vortex
                              // 2 : decaying turbulence
  isgn = +1;
  if (para[5] < 0)
    isgn = -1;

  sigma = para[6]; // initial blob or vortex width

  amp = para[7]; // initial amplitude (of flow, vortex or "turbulent" bath)

  pen = para[8]; // obstacle penalty

  lx = para[9];       // length of y domain in units of rho_s
  nx = int(para[10]); // grid points in x
  if (incon == 0)
    nx *= 2;

  ny = int(para[11]); // grid points in y

  hy = para[10] / lx; // finite-differencing factor, assumes square grid

  wsrc = sigma * double(ny) / lx;

  dt = para[12];         // time step
  itstp = int(para[13]); // inner steps between outputs
  itmax = int(para[14]); // total steps until end

  isolv = int(para[16]);   // Poisson solver
  iadvect = int(para[17]); // advection solver
  dhyv = para[18];         // hyperviscosity coefficient

  npar = int(para[19]); // number of parallel OpenMP threads (depending on available cores)

  // Initialise Open_MP ...................
#ifdef _OPENMP
  ith = fftw_init_threads();
  if (ith == 0)
  {
    printf("| ith=%d : thread init failed... exiting\n", ith);
    exit(1);
  }
  omp_set_num_threads(npar);
  int np, kp;
  printf("| parallel threads ");
#pragma omp parallel private(kp)
  {
    np = omp_get_num_threads();
    kp = omp_get_thread_num();
    printf("%d/", kp);
  }
  printf(" of %d active.\n", np);
  double dnp = double(np);
#else
  printf("| single core processing\n");
  double dnp = 1.;
#pragma omp barrier
#endif

  // Initialise parameters ..................

  nxh = nx / 2;
  nyh = ny / 2;
  nx1 = nx - 1;
  ny1 = ny - 1;

  hy2 = hy * hy;
  hysq = hy2 / 12.;
  xyz = 1. / double(nx * ny);

  if (istart != 2)
  {
    f = fopen("eng.dat", "w");
    fclose(f);
    f = fopen("time.dat", "w");
    fclose(f);
  }

  // coefficient pre-calculation for FFTW Poisson solver
  int ik, jk;
  double cnx = 1. / double(nx), cny = 1. / double(ny);
  double cxy = cnx * cny;
  double kkqq;

  for (int i = 0; i <= nx1; ++i)
    for (int j = 0; j <= ny1; ++j)
    {
      ik = (i >= nxh) ? i - nx : i;
      jk = (j >= nyh) ? j - ny : j;
      // k^2:
      kkqq = (ZwPi * hy * ik * cnx) * (ZwPi * hy * ik * cnx) + (ZwPi * hy * jk * cny) * (ZwPi * hy * jk * cny);
      cvor[i][j] = -cxy / kkqq; // for Poisson equation solver: - k^2 F(k) = Q(k)
    }
  cvor[0][0] = 0.;
}

// ---------------------------------------------------------------------------
// calculate and write output quantities:
// e.g. energies, transport, k spectra, 2D arrays, profiles, ...

// diagnose(ttt, t00, it, itmax, phi, www, www);
void diagnose(double ttt, double t00, int it, int itmax,
              double (*pp)[nymax], double (*ww)[nymax], double (*dd)[nymax], int counter, double reynold, double hypervisc, double dt, double itstp )
{
  int i, j, ik, im, ip, jm, jp;
  char filename_e_mat[150]; // filename buffer always must be defined higher than the _OPENMP stuff
  char filename_w2d[32];
  char filename_v2d[32];
  double enk, enw, dum;
  FILE *f, *g, *h, *g1, *g2, *g3, *g4, *g5, *ene, *grid; // ene is the file where I save the energies at all points

  enk = 0.;
  enw = 0.;
  // energetic quantities
  #ifdef _OPENMP
  #pragma omp parallel private(j, im, ip, jm, jp, dum) shared(enw, enk)
  #endif

  // Save the energy of the gird points into a matrix
  
  //
  for (i = nx0; i <= nx1; i++) //nx0 =nx/2 bcs pipeflow is defined like that here
  {
    // Definiere im und ip sodass differenz 2 zwischen ihnen ist und man so die 
    // Ableitung berechnen kann
    im = (i == nx0) ? nx1 : i - 1;
    ip = (i == nx1) ? nx0 : i + 1;
    for (j = 0; j <= ny1; j++)
    { 
      
      // define a enk_i so u can get energy of each grid point
      jm = (j == 0) ? ny1 : j - 1;
      jp = (j == ny1) ? 0 : j + 1;
      // total enstrophy:
      enw += ww[i][j] * ww[i][j];
      // es gilt  hy = ny / l--> somit Definition der finiten Differenzen
      // total kinetic flow energy:
      vx = .5 * hy * (pp[i][jp] - pp[i][jm]); // vx
      enk += vx * vx; //vx²
      vy = .5 * hy * (pp[ip][j] - pp[im][j]); // vy
      enk += vy * vy;

    }
  }


  grid = fopen("simdata1/info.dat", "w"); // open the grid data file
  fprintf(grid, "%d  %d  %d  %.6e  %.6e  %.6e  %.15e \n", nx/2, ny, counter, reynold, hypervisc,dt,itstp); // bcs if it is pipe flow the grid in x is automatically used by double the value of the inp paramter
  fclose(grid);
  #ifdef _OPENMP
  #pragma omp barrier
  #endif
//ss
  enk *= .5 * xyz;
  enw *= .5 * xyz * xyz;
  if (enw == 0.)
    enw = 1.e-12;
  if (enk == 0.)
    enk = 1.e-12;

  // crash control: stop if nan or inf
  if ((isnan(enw)) || (enw < 1.e-16))
  {
    printf("\n|t=%.3f (%.2f): enn=%.3e\n", ttt, (t00 + itmax * dt), enw);
    printf("|VOR2 CRASHED!  enn=%.3e,   eno=%.3e\n", enw, eno);
  #ifdef _OPENMP
      fftw_cleanup_threads();
  #endif
    exit(1);
  }
  eno = enw;

  // global energy time series output:
  f = fopen("eng.dat", "a");
  if (ttt > 0.)
    fprintf(f, "%.8f  %.6e  %.6e\n", ttt, enw, enk);
  fclose(f);

  // x-cut at y=1:
  g = fopen("cutx0.dat", "w");
  // for (i=0; i<=nx1; i++)
  for (i = nx0; i <= nx1; i++)
  // for (i=nxh; i<=nx1; i++)
  {
    im = (i == 0) ? imbc : i - 1;
    ip = (i == nx1) ? ipbc : i + 1;
    for (j = 1; j <= 1; j++)
    {
      jm = (j == 0) ? 0 : j - 1;
      jp = (j == ny1) ? ny1 : j + 1;
      vx = +.5 * hy * (pp[i][jp] - pp[i][jm]);
      vy = -.5 * hy * (pp[ip][j] - pp[im][j]);
      fprintf(g, "%d  %.6e   %.6e   %.6e\n", i, ww[i][j], pp[i][j], vy);
    }
  }
  fclose(g);

  // x-cut at y=ny/2:
  g = fopen("cutx.dat", "w");
  // for (i=0; i<=nx1; i++)
  for (i = nx0; i <= nx1; i++)
  // for (i=nxh; i<=nx1; i++)
  {
    im = (i == 0) ? imbc : i - 1;
    ip = (i == nx1) ? ipbc : i + 1;
    // im = (i==0) ?   imbc : i-1;
    // ip = (i==nx1) ? ipbc : i+1;
    for (j = nyh; j <= nyh; j++)
    {
      jm = (j == 0) ? 0 : j - 1;
      jp = (j == ny1) ? ny1 : j + 1;
      vx = +.5 * hy * (pp[i][jp] - pp[i][jm]);
      vy = -.5 * hy * (pp[ip][j] - pp[im][j]);
      fprintf(g, "%d  %.6e   %.6e   %.6e\n", i, ww[i][j], pp[i][j], vy);
    }
  }
  fclose(g);

  // x profile (y averages):
  g = fopen("profx.dat", "w");
  snprintf(filename_e_mat,sizeof(char) * 150,"simdata1/e_mat%i.dat",counter);
  ene = fopen(filename_e_mat, "w"); // open the energy file
  double xp_vy;
 
  // for (i=0; i<=nx1; i++)
  for (i = nx0; i <= nx1; i++)
  // for (i=nxh; i<=nx1; i++)
  { 

    im = (i == 0) ? 0 : i - 1;
    ip = (i == nx1) ? nx1 : i + 1;
    xp_vy = 0.;
    for (j = 0; j <= ny1; j++)
    {
      enk_i=0;
      jm = (j == 0) ? 0 : j - 1;
      jp = (j == ny1) ? ny1 : j + 1;
      vx = +.5 * hy * (pp[i][jp] - pp[i][jm]);
      enk_i = vx*vx;
      vy = -.5 * hy * (pp[ip][j] - pp[im][j]);
      enk_i += vy*vy;
      xp_vy += vy;
      // Wichtig fürs verständnis: checke ableitungen bei Zeile ca 840
      fprintf(ene, "%.8e  %d  %d \n", enk_i, i, j); // print the kinetic energie of each grid point into file ene
    }
    fprintf(g, "%d  %.6e\n", i, xp_vy / double(ny));
  }
  fclose(ene); // close the energy file
  fclose(g);

  // y profile (x averages):
  g = fopen("profy.dat", "w");
  for (j = 0; j <= ny1; j++)
  {
    jm = (j == 0) ? 0 : j - 1;
    jp = (j == ny1) ? ny1 : j + 1;
    xp_vy = 0.;
    // for (i=nxh; i<=nx1; i++)
    for (i = nx0; i <= nx1; i++)
    {
      im = (i == 0) ? 0 : i - 1;
      ip = (i == nx1) ? nx1 : i + 1;
      vx = +.5 * hy * (pp[i][jp] - pp[i][jm]);
      vy = -.5 * hy * (pp[ip][j] - pp[im][j]);
      xp_vy += vy;
    }
    fprintf(g, "%d  %.6e\n", j, xp_vy / double(nxh));
  }
  fclose(g);

  // 2D (x,y) plots of vorticity, potential
  snprintf(filename_w2d,sizeof(char) * 32,"simdata1/w2d%i.dat",counter);
  snprintf(filename_v2d,sizeof(char) * 32,"simdata1/v2d%i.dat",counter);
  g1 = fopen(filename_w2d, "w");
  g2 = fopen("p2d0.dat", "w");
  g3 = fopen(filename_v2d, "w");
  int ii;
  double wout, pout, vout, wxx = 0., pxx = 0., vxx = 0.;
  for (i = nxh; i <= nx1; i++)
    for (j = 1; j <= ny1 - 1; j++)
    {
      if (fabs(ww[i][j]) > wxx)
        wxx = fabs(ww[i][j]);
      if (fabs(pp[i][j]) > pxx)
        pxx = fabs(pp[i][j]);
    }

  for (i = nx0; i <= nx1; i++)
  {
    im = (i == 0) ? 0 : i - 1;
    ip = (i == nx1) ? nx1 : i + 1;
    ii = nx1 - i; // x-Achse nach "oben", y-Achse nach "rechts" z.B. in gnuplot
    for (j = 0; j <= ny1; j++)
    {
      jm = (j == 0) ? 0 : j - 1;
      jp = (j == ny1) ? ny1 : j + 1;
      vx = +.5 * hy * (pp[i][jp] - pp[i][jm]);
      vy = -.5 * hy * (pp[ip][j] - pp[im][j]);
      vout = vy - v0avg;
      // vout = sqrt(vx*vx+vy*vy);
      wout = ww[i][j];
      pout = pp[i][j];

      // stamp max/min values for gnuplot colour scale centered to zero
      if (incon == 0)
      {
        if ((i == nxh) && (j == ny1))
        {
          wout = wxx;
          pout = pxx;
        }
        if ((i == nx1) && (j == ny1))
        {
          wout = -wxx;
          pout = -pxx;
        }
      }

      fprintf(g1, "%d  %d  %.6e\n", j, ii, wout);
      fprintf(g2, "%d  %d  %.6e\n", j, ii, pout);
      fprintf(g3, "%d  %d  %.6e\n", j, ii, vout);
    }
    fprintf(g1, "\n");
    fprintf(g2, "\n");
    fprintf(g3, "\n");
  }
  fclose(g1);
  fclose(g2);
  fclose(g3);
  rename("w2d0.dat", "w2d.dat");
  rename("p2d0.dat", "p2d.dat");
  rename("v2d0.dat", "v2d.dat");

  // Fourier ky spectrum
  double py[ny];
  fftw_complex ky[ny / 2], sumky[ny / 2];
  #ifdef _OPENMP
    fftw_plan_with_nthreads(npar);
  #endif
  fftw_plan hindfty;

  hindfty = fftw_plan_dft_r2c_1d(ny, &py[0], &ky[0], FFTW_ESTIMATE);

  for (j = 0; j <= ny / 2; j++)
  {
    sumky[j][0] = 0.;
  }

  for (i = nx0; i <= nx1; i++)
  {
    // for (j=0; j<=ny1; j++) py[j] = ww[i][j];
    // fixed i, iterate over all j
    for (j = 0; j <= ny1; j++)
      py[j] = pp[i][j]; // assign value to py
    // make FT
    fftw_execute(hindfty);
    for (j = 0; j <= ny / 2; j++)
      sumky[j][0] += ky[j][0] * ky[j][0]; // hier werden alle y-FT für verschiedene x summiert
  }
  for (j = 0; j <= ny / 2; j++)
    sumky[j][0] /= (nx * ny);
  fftw_destroy_plan(hindfty);

  g = fopen("pky.dat", "w");
  for (j = 1; j < ny / 2; j++)
    fprintf(g, "%.6e  %.6e\n", hy * double(j) / lx, sumky[j][0]);
  fclose(g);

  // timestamp output
  printf("|t=%.3e (%.3f): enn=%.3e\n", ttt, (t00 + itmax * dt), enw);
//--------------------------------------------------------------------------------------------------------------------------------
  // For Bachelor thesis create a E_k spectrum Christoph Oberladstätter 25.März.2022
  double Esy[nx], sumEsy[nx]; // Energie im Ortsraum (Space)
  fftw_complex Eky[nx/2 ], sumEky[nx /2];
  #ifdef _OPENMP
    fftw_plan_with_nthreads(npar);
  #endif
  fftw_plan hindftEsy;

  hindftEsy = fftw_plan_dft_r2c_1d(nx, &Esy[0], &Eky[0], FFTW_ESTIMATE);
  for (i = 0; i <= nx/2; i++)
  {
    sumEky[i][0] = 0.;
  }
double vxc, vyc;
  for (j = 80; j <= 85; j++) // sum over this y-range
  {
      jm = (j == 0) ? ny1 : j - 1;
      jp = (j == ny1) ? 0 : j + 1;
      for (i = 0; i <= nx1; i++)
      {
      enk=0;
      Esy[i]=0;
      // Definiere im und ip sodass differenz 2 zwischen ihnen ist und man so die 
      // Ableitung berechnen kann
      im = (i == 0) ? nx1 : i - 1;
      ip = (i == nx1) ? 0 : i + 1;
    
      // total enstrophy:
      //enw += ww[i][j] * ww[i][j];

      // es gilt  hy = ny / l--> somit Definition der finiten Differenzen
      // total kinetic flow energy:
      vxc = .5 * hy * (pp[i][jp] - pp[i][jm]); // vx² ----- ist es nicht nur vx??
      vyc = .5 * hy * (pp[ip][j] - pp[im][j]); // vy²
      enk = vxc*vxc + vyc*vyc;
      Esy[i] = enk;
      }
    fftw_execute( hindftEsy ); // Esy--> Eky
    for(i = 0; i <= nx/2; i++)
    {
      sumEky[i][0] += Eky[i][0]; // sum all the Eky to the sumEky
      
    }
  }

  for (i = 0; i <= nx/2; i++)
  {
    sumEky[i][0] /= (nx * ny)*5; // norm the Fourier transfrom---> AAAAAAAAAAAAAAACHTUNG hier wurde 5 für 200-205 hinzugefügt
  }

  fftw_destroy_plan(hindftEsy);  
  g4 = fopen("sumEky.dat", "w");
  g5= fopen("sumEsy.dat", "w");

  for (i = 1; i <= nx/2; i++)
  {
    fprintf(g4, "%.6e  %.6e\n",  double(i) / lx, sqrt(sumEky[i][0]*sumEky[i][0] + sumEky[i][1] * sumEky[i][1]));
    fprintf(g5, "%.6e  %.6e\n", double(i) / lx, sumEsy[i]);
  }
  fclose(g4);
  fclose(g5);
}

// ---------------------------------------------------------------------------
// symmetric mirroring in x for better FFT
void mirror(double (*arrin)[nymax])
{
  int i, j, i0 = 0, ict, jct;
#ifdef _OPENMP
#pragma omp parallel for private(j, ict, jct)
#endif
  for (i = i0; i < nxh + 1; i++)
  {
    ict = nx1 - i + i0;
    for (j = 0; j <= ny1; ++j)
    {
      jct = j;
      arrin[i][j] = arrin[ict][jct];
    }
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

// ---------------------------------------------------------------------------
// Arakawa 4th order scheme for advective 2D Poisson brackets
// [u,y] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
// ideally conserving energy and enstrophy
// [A. Arakawa, Journal of Computational Physics, 119 (1966)]
void arakawa(double (*uuu)[nymax], double (*vvv)[nymax], double (*ww)[nymax])
{
  int i0, j0, ip, jp, im, jm;
  double xxx;

#ifdef _OPENMP
#pragma omp parallel for private(im, ip, j0, jm, jp, xxx)
#endif
  for (i0 = 0; i0 <= nx1; i0++)
  {
    im = (i0 == 0) ? imbc : i0 - 1;
    ip = (i0 == nx1) ? ipbc : i0 + 1;
    for (j0 = 0; j0 <= ny1; ++j0)
    {
      jm = (j0 == 0) ? jmbc : j0 - 1;
      jp = (j0 == ny1) ? jpbc : j0 + 1;
      xxx = vvv[i0][jm] * (+uuu[ip][j0] - uuu[im][j0] - uuu[im][jm] + uuu[ip][jm]);
      xxx += vvv[i0][jp] * (-uuu[ip][j0] + uuu[im][j0] - uuu[ip][jp] + uuu[im][jp]);
      xxx += vvv[ip][j0] * (+uuu[i0][jp] - uuu[i0][jm] + uuu[ip][jp] - uuu[ip][jm]);
      xxx += vvv[im][j0] * (-uuu[i0][jp] + uuu[i0][jm] + uuu[im][jm] - uuu[im][jp]);
      xxx += vvv[ip][jm] * (+uuu[ip][j0] - uuu[i0][jm]);
      xxx += vvv[ip][jp] * (+uuu[i0][jp] - uuu[ip][j0]);
      xxx += vvv[im][jm] * (+uuu[i0][jm] - uuu[im][j0]);
      xxx += vvv[im][jp] * (+uuu[im][j0] - uuu[i0][jp]);
      ww[i0][j0] = hysq * xxx;
    };
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

// ---------------------------------------------------------------------------
// Simple 2nd order upwind scheme for advective 2D Poisson brackets in pipe flow
// [u,y] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
void upwind(double (*uuu)[nymax], double (*vvv)[nymax], double (*ww)[nymax])
{
  int i1, i0, j0, ip, jp, im, jm, im2, ip2, jm2, jp2;
  double xxx, vvx, vvy, ddx, ddy;
  i1 = nx1;
#ifdef _OPENMP
#pragma omp parallel for private(im, ip, im2, ip2, j0, jm, jp, jm2, jp2, xxx)
#endif
  for (i0 = 0; i0 <= i1; i0++)
  {
    im = (i0 == 0) ? imbc : i0 - 1;
    ip = (i0 == i1) ? ipbc : i0 + 1;
    im2 = (im == 0) ? imbc - 1 : im - 1;
    ip2 = (ip == i1) ? ipbc + 1 : ip + 1;

    for (j0 = 0; j0 <= ny1; ++j0)
    {
      jm = (j0 == 0) ? jmbc : j0 - 1;
      jp = (j0 == ny1) ? jpbc : j0 + 1;
      jm2 = (jm == 0) ? jmbc - 1 : jm - 1;
      jp2 = (jp == ny1) ? jpbc + 1 : jp + 1;

      vvx = +.5 * hy * (uuu[i0][jp] - uuu[i0][jm]);
      vvy = -.5 * hy * (uuu[ip][j0] - uuu[im][j0]);

      // for test: only 1st order upwind in y pipe, but centered in x:
      // xxx = .5*(uuu[ip][j0]-uuu[im][j0])*(vvv[i0][j0]-vvv[i0][jm]);
      // xxx-= .5*(vvv[ip][j0]-vvv[im][j0])*(uuu[i0][j0]-uuu[i0][jm]);

      // 2nd order upwind
      if (vvx > 0)
      {
        ddx = 3. * uuu[i0][j0] - 4. * uuu[im][j0] + uuu[im2][j0];
      }
      else
      {
        ddx = 4. * uuu[ip][j0] - 3. * uuu[i0][j0] - uuu[ip2][j0];
      };
      if (vvy > 0)
      {
        ddy = 3. * vvv[i0][j0] - 4. * vvv[i0][jm] + vvv[i0][jm2];
      }
      else
      {
        ddy = 4. * vvv[i0][jp] - 3. * vvv[i0][j0] - vvv[i0][jp2];
      };
      xxx = .25 * ddx * ddy;
      if (vvx > 0)
      {
        ddx = 3. * vvv[i0][j0] - 4. * vvv[im][j0] + vvv[im2][j0];
      }
      else
      {
        ddx = 4. * vvv[ip][j0] - 3. * vvv[i0][j0] - vvv[ip2][j0];
      };
      if (vvy > 0)
      {
        ddy = 3. * uuu[i0][j0] - 4. * uuu[i0][jm] + uuu[i0][jm2];
      }
      else
      {
        ddy = 4. * uuu[i0][jp] - 3. * uuu[i0][j0] - uuu[i0][jp2];
      };
      xxx -= .25 * ddx * ddy;

      ww[i0][j0] = -hy2 * xxx;
    };
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

// ---------------------------------------------------------------------------
// Simple 2nd order centered difference scheme for advective 2D Poisson bracket
// [u,y] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
void simple(double (*uuu)[nymax], double (*vvv)[nymax], double (*ww)[nymax])
{
  int i0, j0, ip, jp, im, jm;
  double xxx;
#ifdef _OPENMP
#pragma omp parallel for private(im, ip, j0, jm, jp, xxx)
#endif
  for (i0 = 0; i0 <= nx1; i0++)
  {
    im = (i0 == 0) ? imbc : i0 - 1;
    ip = (i0 == nx1) ? ipbc : i0 + 1;
    for (j0 = 0; j0 <= ny1; ++j0)
    {
      jm = (j0 == 0) ? jmbc : j0 - 1;
      jp = (j0 == ny1) ? jpbc : j0 + 1;
      xxx = (uuu[ip][j0] - uuu[im][j0]) * (vvv[i0][jp] - vvv[i0][jm]);
      xxx -= (vvv[ip][j0] - vvv[im][j0]) * (uuu[i0][jp] - uuu[i0][jm]);
      ww[i0][j0] = -hy2 * .25 * xxx;
    };
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

// ---------------------------------------------------------------------------
// 2D Laplace operator by finite differences: y = del^2 x
// standard: 2nd order 9pt stencil
// (for 4th order version: uncomment the respective code lines)
void laplace(double (*fo)[nymax], double (*fi)[nymax], int ib, int io)
{
  int i, j, im, ip, jm, jp, i0, i1;
  int im2, ip2, jm2, jp2;

  i0 = nx0 + ib;
  i1 = nx1 - ib;
  i0 = 0;
  i1 = nx1;

  if (ib > 0)
    i0 = nxh + ib;

#ifdef _OPENMP
#pragma omp parallel for private(im, ip, im2, ip2, j, jm, jp, jm2, jp2)
#endif
  for (i = i0; i <= i1; i++)
  {
    im = (i == i0) ? imbc : i - 1;
    ip = (i == i1) ? ipbc : i + 1;

    // im2 = (im==i0) ? i1-1 : im-1; // for higher order / periodic
    // ip2 = (ip==i1) ? i0+1 : ip+1; // for higher order / periodic

    for (j = 0; j <= ny1; j++)
    {
      jm = (j == 0) ? jmbc : j - 1;   // for in/out-flow b.c.
      jp = (j == ny1) ? jpbc : j + 1; // for in/out-flow b.c.

      // jm2 = (jm==0)   ? ny1-1 : jm-1; // for higher order / periodic
      // jp2 = (jp==ny1) ? 0+1   : jp+1; // for higher order / periodic

      // 5-pt stencil (2nd order)
      /*
      fo[i][j] = fi[ip][j]+fi[im][j]+fi[i][jp]+fi[i][jm]-4.*fi[i][j];
      fo[i][j]*= hy2;
      // */

      // 9-pt stencil (2nd order)
      // /*
      fo[i][j] = 4. * fi[ip][j] + 4. * fi[im][j] + 4. * fi[i][jp] + 4. * fi[i][jm] + fi[im][jm] + fi[im][jp] + fi[ip][jm] + fi[ip][jp] - 20. * fi[i][j];
      fo[i][j] *= hy2 / 6.;
      // */

      // 4th order version requiring two neighbours on each side:
      /*
      fo[i][j] = - 60.*fi[i][j] + 16.*(fi[ip][j]+fi[i][jp]+fi[im][j]+fi[i][jm])
        - (fi[ip2][j]+fi[im2][j]+fi[i][jp2]+fi[i][jm2]);
      fo[i][j]*= hysq;
      // */
    }
  }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

// ---------------------------------------------------------------------------
// 2nd order SOR iterative 2D Poisson equation solver
// equation: " laplace(uuu) = qqq " : gets uuu  for given qqq
void sor2p(double (*qqq)[nymax], double (*uuu)[nymax])
{
  int i, j, k, im, ip, jm, jp, i0, i1, j0, j1, iter, nit, idx, ibb, lll, dum, i0p, i1m;
  double dx, dy, dxy, dxy2;
  double resid, rel, rel2, pemit, epsmax, om_sor, fcheck, epssum, fnorm = 0.;
  double uu1p[nx][ny], fffp[nx][ny];

  i0 = 0;
  i1 = nx1;
  j0 = 0;
  j1 = ny1;

  nit = 500;      // max. number of iterations
  epsmax = 0.001; // max. allowed error

  ibb = 0;
  i0p = i0 + ibb;
  i1m = i1 - ibb;

  idx = i1 - i0;
  rel = (cos(2. * M_PI / double(idx)) + cos(1. * M_PI / double(ny))) / 2.; // mix b.c.
  // rel = ( cos(2.*M_PI/double(idx)) + cos(2.*M_PI/double(ny)) ) / 2.; // per.
  rel2 = .25 * rel * rel;

  // preparation of coefficients
#ifdef _OPENMP
#pragma omp parallel for private(im, ip, j, jm, jp)
#endif
  for (i = i0; i <= i1; ++i)
  {
    im = (i == i0) ? imbc : i - 1;
    ip = (i == i1) ? ipbc : i + 1;
    for (j = j0; j <= j1; ++j)
    {
      jm = (j == j0) ? jmbc : j - 1;
      jp = (j == j1) ? jpbc : j + 1;

      uu1p[i][j] = 0.; // initial value for iteration
      if (isorcnt > 0)
        uu1p[i][j] = uuu[i][j];
      if (isorcnt > 1)
        uu1p[i][j] += (uuu[i][j] - uu0[i][j]);

      fffp[i][j] = qqq[i][j];
      // 2nd order 5pt correction to r.h.s.
      // -> 4th order discretization of Laplacian (else only 2nd order)
      fffp[i][j] += (qqq[im][j] + qqq[ip][j] + qqq[i][jm] + qqq[i][jp] - 4. * qqq[i][j]) * (1. / 12.);
      fnorm += fabs(fffp[i][j]);
      fffp[i][j] *= 1. / hy2;
    }
  }
#ifdef _OPENMP
#pragma omp barrier
#endif

  isorcnt++;
  fcheck = epsmax * fnorm;
  om_sor = 1.;

  // begin iteration
  for (iter = 1; iter <= nit; ++iter)
  {
    epssum = 0.;
    // even/odd ordering
    for (lll = 0; lll <= 1; lll++)
    {
      dum = lll;
#ifdef _OPENMP
#pragma omp parallel for private(im, ip, j, jm, jp, resid)
#endif
      for (i = i0p; i <= i1m; i++) // ibb=1 -> phi=0 b.c.
      {
        im = (i == i0) ? imbc : i - 1;
        ip = (i == i1) ? ipbc : i + 1;
        for (j = j0 + dum; j <= j1; j += 2)
        {
          jm = (j == j0) ? jmbc : j - 1;
          jp = (j == j1) ? jpbc : j + 1;

          resid = (uu1p[ip][j] + uu1p[im][j] + uu1p[i][jp] + uu1p[i][jm] - 4. * uu1p[i][j] - fffp[i][j]);
          uu1p[i][j] += .25 * om_sor * resid;
          epssum += fabs(resid);
        }
        dum = 1 - dum;
      }
#ifdef _OPENMP
#pragma omp barrier
#endif
      om_sor = ((iter == 1) && (lll == 0)) ? 1. / (1. - 2. * rel2) : 1. / (1. - rel2 * om_sor);
    }
    if (epssum < fcheck)
    {
      break;
    }
  }
  // end iteration

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = i0; i <= i1; ++i)
    for (j = j0; j <= j1; ++j)
    {
      uu0[i][j] = uuu[i][j];
      uuu[i][j] = uu1p[i][j];
    }
#ifdef _OPENMP
#pragma omp barrier
#endif
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
