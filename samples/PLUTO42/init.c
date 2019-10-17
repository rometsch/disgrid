/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet interaction problem.

  Simulate the interaction of a planet embedded in a disk as described
  in section 3.4 of Mignone et al., A&A (2012) 545, A152.
  This test is a nice benchmark for the FARGO module and the
  \c ROTATING_FRAME switch.
  For testing-purposes  no viscosity is used here.
  The initial condition consists of a locally isothermal configuration
  with temperature profile \f$\propto T^{-1}\f$  yielding a disk vertical
  height to radius of \c 0.05.
  The gravitational potential due to the presence of the star and the
  planet is defined in BodyForcePotential() function.

  The conventions used throught the implementation are the following:
 
  - \c r  = spherical radius
  - \c R  = cylindrical radius
  - \c z  = cylindrical height
  - \c th = meridional angle
  
  The test can be carried out in polar (2D or 3D) or spherical (3D)
  coordinates and the following parameters determine the initial configuration:
   
  -# <tt>g_inputParam[Mstar]</tt>: controls the star mass (in solar masses)
  -# <tt>g_inputParam[Mdisk]</tt>: controls the disk mass (in solar masses)
  -# <tt>g_inputParam[Mplanet]</tt>: sets the planet mass (in earth masses)
  -# <tt>g_inputParam[Viscosity]</tt>: sets the amount of viscosity


  Computation can be carried in the rotating or in the observer's frame
  of reference (\c ROTATING_FRAME to \c YES or \c NO, respectively).
  In particular:

  - Configurations #01 and #02 are in 2D polar coordinates without and with
    FARGO, in the rotating frame.
  - Configurations #03, #04 and #05 are in spherical 3D coordinates with
    and without FARGO in the rotating frame
  - Configurations #06 and #07 are in 2D polar coordinates but in the
    observer's frame
  - Configuration #08 employs static AMR  (grid levels are spatial
    dependent but not dependent on time) in the rotating frame.
    

  \image html hd_disk_planet.08.png "Density map for configuration #08 using AMR" width=1cm

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012

  \b References:
     - "A Conservative orbital advection scheme for simulations
        of magnetized shear flows with the PLUTO Code"
        Mignone et al, A&A (2012)
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#ifdef PARALLEL
  #define IF_ROOT if(prank == 0)
#else
  #define IF_ROOT if(1)
#endif

#define UNIT_MASS (UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
#define GM 1.0
//defined in GNU extension of math.h:
void sincos(double x, double *sin, double *cos);

#define MIN_DENSITY 1e-8

static void NormalizeDensity (const Data *d, Grid *g);
#if ROTATING_FRAME == NO
 #define g_OmegaZ  0.0
#endif

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, cs;
  double scrh;

  #if EOS == IDEAL
   g_gamma = 1.01;
  #endif

  #if ROTATING_FRAME == YES
   g_OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
   g_OmegaZ *= 2.0*CONST_PI;
  #endif
  
  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
  #elif GEOMETRY == SPHERICAL
   r  = x1;
   th = x2;
   R  = r*sin(th);
   z  = r*cos(th);
  #endif
  
  H      = 0.05*R;
  OmegaK = 2.0*CONST_PI/(R*sqrt(R));
  cs     = H*OmegaK;
  
  scrh   = (0.5*CONST_PI - th)*r/H;
  us[RHO] = 1.0/(R*sqrt(R))*exp(-0.5*scrh*scrh);
  us[VX1] = us[VX2] = us[VX3] = 0.0;

  us[iVPHI] = R*(OmegaK - g_OmegaZ);
  #if EOS == IDEAL
   us[PRS] = us[RHO]*cs*cs;
  #elif EOS == ISOTHERMAL
//   g_isoSoundSpeed = cs;
   g_isoSoundSpeed = CONST_PI*0.1;
  #endif

#if DUST == YES
  us[RHO_D] = 1.e-4;
  us[VX1_D] = us[VX2_D] = us[VX3_D] = 0.0;
  us[VX2_D] = us[iVPHI];
#endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
{
  double h = 0.05;
  double Hsq,eps2,dV,dx,dy,x,y,xp,yp,d1,d2,d3,dm,dT;
  double sin_phi,cos_phi;
  double mass = 0.0;
  double eng = 0.0;
  double torque = 0.0;
  double Mp = g_inputParam[Mplanet]*CONST_Msun/UNIT_MASS;
  double r,phi,RHill2,sin2;
  double t = g_time;
  double totalmass,totalenergy,totaltorque;
  FILE *out;
  int i,j,k;
  char flnm[] = "out/info.dat";
  double rp = 1.0;

  #if ROTATING_FRAME == YES
    xp = rp;
    yp = 0.0;
  #else
    double phi_p = sqrt(1/(rp*rp*rp))*t;
    xp = rp*cos(phi_p);
    yp = rp*sin(phi_p);
  #endif

  IF_ROOT
  { 
    if (g_time == 0.0)
    {
      out = fopen(flnm,"w");
      fprintf(out, "time\tmass\tenergy\ttorque (excl. 0.5 RHill) \n");
    }
    else out = fopen(flnm,"a");
  }
  
    RHill2 = pow(Mp/3.0, 2.0/3.0);

  KDOM_LOOP(k)
  JDOM_LOOP(j)
  {
    phi = grid[JDIR].x[j];
    IDOM_LOOP(i)
    {
      r = grid[IDIR].x[i];
      dV = grid[IDIR].dV[i] * grid[JDIR].dV[j] * grid[KDIR].dV[k];
      dm = d->Vc[RHO][k][j][i] *dV;

      sincos(phi,&sin_phi,&cos_phi);

      dx = r*cos_phi - xp;      
      dy = r*sin_phi - yp;

      // get smoothing (directly squared)
      #if EOS == IDEAL
        Hsq = r*r*r* d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
      #else
        Hsq = h*h;
      #endif
      eps2 = 0.36* Hsq;

      d2 = dx*dx + dy*dy + eps2;

      if ( d2 > 0.25*RHill2 )
      {
        d1 = sqrt(d2);
        d3 = d1*d2;
        dT = Mp * dm * (xp*dy - yp*dx) /d3;
        torque += dT;
      }

      mass += dm;
      #if EOS != ISOTHERMAL
        eng += (d->Vc[PRS][k][j][i])*dV;
      #endif
    }
  }

  #if EOS != ISOTHERMAL
    eng /= (g_gamma-1.0);
  #endif

  #ifdef PARALLEL
    MPI_Reduce(&mass, &totalmass,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&eng, &totalenergy,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&torque, &totaltorque,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);  //is it really needed?
  #endif

  IF_ROOT
  {
    fprintf(out,"%.16le\t%.16le\t%.16le\t%.16le\n",t,totalmass,totalenergy,totaltorque);
    fflush(out);
    fclose(out);
  }
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3, R, OmegaK, v[256];
  static int do_once = 1;
  
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  #if DIMENSIONS == 3
  if (side == 0){
    if (do_once){
      NormalizeDensity(d, grid);
      do_once = 0;
    }
  }
  #endif

  if (side == X1_BEG){
    X1_BEG_LOOP(k,j,i){
      NVAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] *= -1.0;
      #if GEOMETRY == POLAR
       R = x1[i];
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
#if DUST == YES      
//      NDUST_LOOP(nv) d->Vc[nv][k][j][i] = 0.0;
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
#endif      
    }
  }

  if (side == X1_END){
    X1_END_LOOP(k,j,i){
      NVAR_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];
      #if GEOMETRY == POLAR
       R = x1[i];
//       d->Vc[iVR][k][j][i] = 0.0;
      #elif GEOMETRY == SPHERICAL
       R = x1[i]*sin(x2[j]);
       d->Vc[iVR][k][j][i]  = 0.0;
       d->Vc[iVTH][k][j][i] = 0.0;
      #endif
      OmegaK = 2.0*CONST_PI/(R*sqrt(R));
      d->Vc[iVPHI][k][j][i] = R*(OmegaK - g_OmegaZ);
#if DUST == YES      
      d->Vc[VX2_D][k][j][i] = d->Vc[iVPHI][k][j][i];
#endif      
 
    }
  }
}

/* ************************************************************** */
void NormalizeDensity (const Data *d, Grid *grid)
/*
 *
 * Normalize density and pressure as   rho -> K*rho, where
 *
 *   K = M/(\sum rho*dV)
 *
 **************************************************************** */
{
  int   i, j, k;
  double mc;
        
  mc  = 0.5*g_inputParam[Mdisk]*CONST_Msun;
  mc /= UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH;
  DOM_LOOP(k,j,i){
    d->Vc[RHO][k][j][i] *= mc;
    #if EOS == IDEAL
     d->Vc[PRS][k][j][i] *= mc;
    #endif
  }
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi;

  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
   x  = R*cos(x2);
   y  = R*sin(x2);
  #elif (GEOMETRY == SPHERICAL)
   r  = x1;
   th = x2;
   R = r*sin(th);
   z = r*cos(th);
   x = r*sin(th)*cos(x3);
   y = r*sin(th)*sin(x3);
  #endif

/* ---------------------------------------------
             planet position
   --------------------------------------------- */

  #if ROTATING_FRAME == NO
   double OmegaZ;
   t = g_time;
   if (g_stepNumber == 2) t += g_dt;
   OmegaZ  = sqrt(1.0 + g_inputParam[Mplanet]/g_inputParam[Mstar]*CONST_Mearth/CONST_Msun);
   OmegaZ *= 2.0*CONST_PI;

   xp = cos(OmegaZ*t);
   yp = sin(OmegaZ*t);
  #else
   xp = 1.0/sqrt(2.0);  /* initial planet position */
   yp = 1.0/sqrt(2.0); 
  #endif

  d = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z);
  rsm = 0.03*R;
  if (d > rsm) phiplanet = g_inputParam[Mplanet]/d;
  else phiplanet = g_inputParam[Mplanet]/d*(pow(d/rsm,4.)-2.*pow(d/rsm,3.)+2.*d/rsm);
  
  phi  = - 4.0*CONST_PI*CONST_PI/g_inputParam[Mstar];
  phi *= (g_inputParam[Mstar]/r + phiplanet*CONST_Mearth/CONST_Msun);

  return phi;
}
#endif

