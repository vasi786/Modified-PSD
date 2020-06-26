#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gettime.c"
#include "solvopt.h"
#include "ran_nr.cpp"

#define big_no 1000.0
#define PI 3.141592654
#define MAX_TYPE 120

using namespace std;

long nearest(double x)
{
  if(x)
    {
      return rint(x);
    }
  else
    return 0;
}


class Miscellaneous
{
  public:
    long int RAN_IN;
    Miscellaneous()
    {
      RAN_IN=-3635;
      return;
    }

    double random(double l_limit, double u_limit)
    {
       double core, ret_val;
       //core=genrand();
	   core=ran0(&RAN_IN);
       ret_val=l_limit+core*(u_limit-l_limit);
       return ret_val;
    }

    long int randomInt(long int l_limit, long int u_limit)
    {
       double core, ret_val;
       long int ret_val_int;
       //core=genrand();
	   core=ran0(&RAN_IN);

	   ret_val_int=(long int)(core*(u_limit-l_limit+1))+l_limit;
       if(ret_val_int<l_limit)
          ++ret_val_int;
       if(ret_val_int>u_limit)
          --ret_val_int;
       return ret_val_int;

    }

	double max(double a, double b)
	{
	  if(a>b)
	    return a;
	  else
	    return b;
    }

	double min(double a, double b)
	{
	  if(a<b)
	    return a;
	  else
	    return b;
    }
};



class dLattice
{
public:
  double x;
  double y;
  double z;
  char mol_name[3];
  int mol_type;
  double mol_rad;

  dLattice()
    {
      return;
    }

  dLattice(double x1, double y1, double z1)
    {
      x=x1; y=y1; z=z1;
      return;
    }

  void locate(double x1, double y1, double z1)
    {
      x=x1; y=y1; z=z1;
      return;
    }
};


double lx, ly, lz;
dLattice *solv, *molecules;
long int solvn=0, nmol=0;
double RES_PSD, adsb_rad;
int npoints=solvn;
double xp, yp, zp;
double window=8000;


class psd
{
public:
  double *PSD_hist;
  double *diff_hist;
  double *diff_hist_old;
  double *diff_pt;
  long hist_count;

  psd()
    {
      Miscellaneous misc;
      hist_count=nearest(misc.max(misc.max(lx,ly),lz)*RES_PSD)+1;
      PSD_hist=new double[hist_count];
      diff_hist=new double[hist_count];
      diff_hist_old=new double[hist_count];
      diff_pt=new double[hist_count];
      for(int i=0; i<hist_count; ++i)
	{
	  PSD_hist[i]=0;
	  diff_hist[i]=0;
	  diff_hist_old[i]=0;
	  diff_pt[i]=0;
	}

	 molecules=new dLattice[solvn];
    }

  void find_psd(void);
};

void psd :: find_psd()
{
  double solvopt(unsigned short n,
                 double x[],
                 double fun(double []),
                 double grad(double [], double []),
                 double options[],
                 double func(double []),
                 double gradc(double [], double [])
              );
  Miscellaneous rand;
  npoints=10;
  double *radius=new double[npoints];
  dLattice *points=new dLattice[npoints];
  double fun(double []), func(double []);
  double grad(double [], double []), gradc(double [], double []);
  void build_mol_list(double []);
  double param[3]={0, 0, 0};
  int n=3;
  double f;
  double term=1000;
  long step=1;
  FILE *f_error=fopen("error.dat", "w");
  FILE *f_solution=fopen("solution.log", "w");   //Writing to a file solution.log for further analyses
  build_mol_list(param);

double sec;
double sec1=0;
gettime(&sec);

cout << "Modified by vasista and sandeep" << endl;

while(term)
{
   

  for(int i=0; i<npoints; ++i)
  {
	//select random point
		double dis=0;
		double px=rand.random(0,lx), py=rand.random(0,ly), pz=rand.random(0,lz);
		param[0]=px; param[1]=py; param[2]=pz;
		xp=px; yp=py; zp=pz;
		dis=-fun(param);

	solvopt_options[4]=-1;
  
	if(dis>0)      //not inside one of the spheres and n = 3 dimensions of current system
	  {
      //cout << "random point is " <<param[0] << " " << param[1] << " " << param[2] << endl;
	    f=solvopt(n, param, &fun, (typeof(&grad))(null_entry), solvopt_options, &func, (typeof(&gradc))(null_entry));
      points[i].x=param[0]; points[i].y=param[1]; points[i].z=param[2];
	    radius[i]=sqrt(-f);
      //cout << "radius is " << radius[i] << endl;
       fprintf(f_solution, "%lg %lg %lg %lg\n", param[0], param[1],param[2],radius[i]);
       fflush(f_solution);
	  }
    
	else
	  {
	    radius[i]=-1;
	  }
  }
    



  //gettime(&sec1);
  //fprintf(f_stat, "At %lg sec, Creating histogram.....\n", (sec1)-(sec)); fflush(f_stat);
 
  cout << step << "    Writing\n";
  double max_cumm=0; 
  for(int i=0; i<npoints; ++i)
  {
    if(radius[i]>adsb_rad)
      {
		  for(int p=0; p<radius[i]*RES_PSD+1; ++p)
		    {
		      PSD_hist[p]+=1;
		      if(PSD_hist[p]>max_cumm) max_cumm=PSD_hist[p];
		      if(p>=hist_count) cout << p << "\n";
		    }
      }
  }
  //gettime(&sec1);
  //fprintf(f_stat, "At %lg sec, Writing files......\n", (sec1)-(sec)); fflush(f_stat);
 

  double max_diff=0;
  int diff_lim=0;
  diff_pt[0]=0;
  diff_hist[0]=(PSD_hist[0]-PSD_hist[1]);
  max_diff=diff_hist[0];
  diff_pt[1]=1/RES_PSD;
  diff_hist[1]=(PSD_hist[0]-PSD_hist[2])/2;
  if(diff_hist[1]>max_diff) max_diff=diff_hist[1];
  for(int i=2; i<hist_count-2; ++i)
    {
       diff_pt[i]=i/RES_PSD;
       //diff_hist[i]=-(PSD_hist[i-2]-8*PSD_hist[i-1]+8*PSD_hist[i+1]-PSD_hist[i+2])/14;
       diff_hist[i]=(PSD_hist[i-1]-PSD_hist[i+1])/2.0;
       if(diff_hist[i]>max_diff) max_diff=diff_hist[i];
       if(PSD_hist[i]) diff_lim=i+5;
    }

  double err=0, avg_err=0, max_err=0;
  int count=0;
  for(int i=0; i<hist_count-2; ++i)
    {
      diff_hist[i]/=max_diff;
      if(diff_hist[i])
	{
	  err=fabs(diff_hist[i]-diff_hist_old[i])/diff_hist[i];
	  avg_err+=err;
	  if(err>max_err) max_err=err;
	  ++count;
	}
    }
  avg_err/=count;
  gettime(&sec1);
  fprintf(f_error, "%ld %lg %lg %lg\n", step, avg_err, max_err, (sec1-sec)/3600);
  fflush(f_error);

  FILE *f_psd_cumm=fopen("psd_cumm.dat", "w");
  for(int i=0; i<diff_lim; ++i)
    {
      fprintf(f_psd_cumm, "%lg %lg\n", i/RES_PSD, PSD_hist[i]/max_cumm);
    }
  fclose(f_psd_cumm);
  
  FILE *f_psd_diff=fopen("psd_diff.dat", "w");
  for(int i=0; i<hist_count-2; ++i)
    {
      if(i<=diff_lim) fprintf(f_psd_diff, "%lg %lg\n", diff_pt[i], diff_hist[i]);
      diff_hist_old[i]=diff_hist[i];
    }
  fclose(f_psd_diff);

  ++step;

}
};



main(int argc, char **argv)
{
  if(argc<7)
  {
	cout << "Usage: psd < xyz-file > < test particle radius(A) > < bin-size(A) > < lx > < ly > < lz >\n";
	exit(-1);
  }
  FILE *f_cris=fopen(argv[1], "r");
  adsb_rad=strtod(argv[2], NULL);
  RES_PSD=1/strtod(argv[3], NULL);
  lx=strtod(argv[4],NULL);
  ly=strtod(argv[5],NULL);
  lz=strtod(argv[6],NULL);
  char dummy, header[100];
  long int datan;
  fscanf(f_cris, "%ld", &datan);
  fscanf(f_cris, "%c%[^\n]", &dummy, header);
  dLattice *data=new dLattice[datan];
  solv=new dLattice [datan];
  solvn=0;
  for(long int i=0; i<datan; ++i)
    {
      fscanf(f_cris, "%c%s %lg %lg %lg", &dummy, data[i].mol_name, &data[i].x, &data[i].y, &data[i].z);
      //if(data[i].mol_name=='S' || data[i].mol_name=='O')  // || data[i].mol_name=='H')
	{
	  solv[solvn]=data[i];
	  ++solvn;
	}
    }
  fclose(f_cris);
  delete[] data;

  dLattice box_min=dLattice(1e6, 1e6, 1e6), box_max=dLattice(-1e6, -1e6, -1e6);
  for(int i=0; i<solvn; ++i)
    {
      if(solv[i].x<box_min.x) box_min.x=solv[i].x;
      if(solv[i].y<box_min.y) box_min.y=solv[i].y;
      if(solv[i].z<box_min.z) box_min.z=solv[i].z;
      if(solv[i].x>box_max.x) box_max.x=solv[i].x;
      if(solv[i].y>box_max.y) box_max.y=solv[i].y;
      if(solv[i].z>box_max.z) box_max.z=solv[i].z;
    }

  //lx=box_max.x-box_min.x; ly=box_max.y-box_min.y; lz=box_max.z-box_min.z;
    //lx=68.565002;
    //ly=64.458000;
    //lz=81.552002;
  for(int i=0; i<solvn; ++i)
  {
	  solv[i].x-=box_min.x;
	  solv[i].y-=box_min.y;
	  solv[i].z-=box_min.z;
  }
  printf("Box Dimensions (A): %lg x %lg x %lg\n", lx, ly, lz);
  //Read radii list
  FILE *f_rad=fopen("radii_list.dat", "r");
  if(!f_rad) { cout << "You need a file called radii_list.dat with all the atom types and vdW radii.\nExiting\n"; exit(-1);}
  
  int ntype=0;
  struct
  {
	char type[3];
	double rad;
  }type_list[MAX_TYPE];
  
  while(getc(f_rad)!='\n')
  {
	fseek(f_rad, -1, SEEK_CUR);
	fscanf(f_rad, "%s%lg%[^\n]", type_list[ntype].type, &type_list[ntype].rad, header);
	getc(f_rad);
	++ntype;
  }
  fclose(f_rad);
  cout << "Atomtyping....";
  int found=0;
  for(int i=0; i<solvn; ++i)
  {
	found=0;
	for(int j=0; j<ntype && !found; ++j)
	{
		if(strcmp(solv[i].mol_name, type_list[j].type) == 0)
		{
			solv[i].mol_rad=type_list[j].rad;
			found=1;
		}
	}
	if(!found)
	{
		cout << "Unknown atom type encountered.\nPlease add the atom symbol and radius in res_list.dat\n";
		exit(-1);
	}
  }
  cout << "done\n";
  psd *PSD=new psd();
  PSD->find_psd();
  delete PSD;
  delete[] solv;
}


void build_mol_list(double param[])
{
	int count=0;
	for(int i=0; i<solvn; ++i)
	{

		double x1=solv[i].x, y1=solv[i].y, z1=solv[i].z;	
	   
		double rx=fabs(param[0]-x1), ry=fabs(param[1]-y1), rz=fabs(param[2]-z1);
		if(rx>0.5*lx) rx-=lx;
		if(ry>0.5*ly) ry-=ly;
		if(rz>0.5*lz) rz-=lz;
		double dis=(rx*rx+ry*ry+rz*rz);
		if(dis<=window*window)
		{
			molecules[count]=solv[i];
			++count;
		}
	}
	nmol=count;
}


double fun( double param[])
{
	double minr2=1e10;
	for(int i=0; i<nmol; ++i)
	{
		double x1=molecules[i].x, y1=molecules[i].y, z1=molecules[i].z;
		
		double rx=fabs(param[0]-x1), ry=fabs(param[1]-y1), rz=fabs(param[2]-z1);
		//periodic boundary convention
		if(rx>0.5*lx) rx-=lx;
		if(ry>0.5*ly) ry-=ly;
		if(rz>0.5*lz) rz-=lz;
		double dis=sqrt(rx*rx+ry*ry+rz*rz)-solv[i].mol_rad;
		if(dis<minr2) minr2=dis;
	}
	return minr2>0?-minr2*minr2:0;
}


double func(double param[])
{
        Miscellaneous misc=Miscellaneous();
	double minr2=-fun(param);
	double maxcon=-1e10;
	
	double dis=(xp-param[0])*(xp-param[0])+(yp-param[1])*(yp-param[1])+(zp-param[2])*(zp-param[2]);
	  maxcon=-minr2+dis;

	  return misc.max(0,maxcon);
}













