#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <complex>
#include "omp.h"
#include "mkl.h"
#include "/data1/home/stephan/Photoionization_SAE_PW/version_3.0/photoionization_module_1.0/pice_module.cpp"

double total_cs(int n_theta,int n_phi,double *theta,double *phi,double *mfpad);
bool compute_mfpad(double x_comp,double y_comp,double z_comp,int n_points,double *mfpad,int kp,std::string dipole_address);
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address);
bool angularly_resolved_dipole_reader(double *prefactor,double *theta,double *phi, double * redipole_x,double * imdipole_x,double * redipole_y,double * imdipole_y,double * redipole_z,double * imdipole_z,int n_points,std::string dipole_address);

int main(int argc, char *argv [])
{
   //THIS CODE COMPUTES THE TDMFPAD USING THE SUDDEN IONIZATION APPROXIMATION FROM A TD WAVE FUNCTION IN THE NEUTRAL MOLECULES. IT USES THE PWAPIC METHOD TO GENERATE THE COUPLING ELEMENTS AND COMPUTES THE MFPAD FOR ALL GEOMETRIES AND ELECTRONIC STATES AT A GIVEN KINETIC ENERGY. IT THE MULTIPLIES THE NUCLEAR WAVE FUNCTIONS BY THE CORRESPONDING MFPADS AND COMPUTES THE TDMFPAD BY COMPUTING THE TOTAL DENSITY OF THE TWO WAVE FUNCTIONS, INTEGRATED OVER NUCLEAR COORDINATE.
   //
   omp_set_num_threads(32);
   using namespace std;
   const double Pi(acos(-1));
   double temp(0);
   double up(0);
   double down(0);
   double rectemp(0);
   double imctemp(0);


   //LOCATION OF THE INPUT FILES
   string state_1_wf_address("/data1/home/stephan/wavepack_results_dir/wvpck_astridpulse_0.02/neut_wf_ceppi_state_3.wvpck");
   string state_2_wf_address("/data1/home/stephan/wavepack_results_dir/wvpck_astridpulse_0.02/neut_wf_ceppi_state_4.wvpck");
   string ionization_coupling_file("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_pice_data_newnorm.h5");//LiH_PICE_R_i_j.txt
   string tdmfpad_address("/data1/home/stephan/wavepack_results_dir/wvpck_astridpulse_0.02/tdmfpad_sudden_st_3_4.txt");
   string anisotropy_address("/data1/home/stephan/wavepack_results_dir/wvpck_astridpulse_0.02/anisotropy_sudden.txt");

   ifstream input;
   ofstream output;
   ofstream output2;

   //DIMENSION OF THE R GRID AND OF THE CONTINUUM
   int grid_size(512);
   int n_theta(20);
   int n_phi(40);
   double kp(sqrt(2.*2./27.211));
   int n_points(n_theta*n_phi);
   int n_times(2000);
   int state_index_1(3);
   int state_index_2(4);
   int state_index_cat(0);

   //ORIENTATION OF THE ELECTRIC FIELD IN THE MOLECULAR FRAME IN SPHERICAL COORD AND TRANSFO TO CARTESIAN
   double theta_elec(0);//(acos(-1)/2);
   double phi_elec(0);//(acos(-1)/2);
   double x_comp(sin(theta_elec)*cos(phi_elec));
   double y_comp(sin(theta_elec)*sin(phi_elec));
   double z_comp(cos(theta_elec));

   //INITIALIZATION OF ALL ARRAYS.
   double *thet=new double[n_points];
   double *phi=new double[n_points];
   double overlap(0);
   double *norm1=new double[n_times];
   double *norm2=new double[n_times];

   std::complex<double> **dipole_1_x=new std::complex<double>*[n_points];
   std::complex<double> **dipole_1_y=new std::complex<double>*[n_points];
   std::complex<double> **dipole_1_z=new std::complex<double>*[n_points];
   std::complex<double> **dipole_2_x=new std::complex<double>*[n_points];
   std::complex<double> **dipole_2_y=new std::complex<double>*[n_points];
   std::complex<double> **dipole_2_z=new std::complex<double>*[n_points];
   double **mfpad=new double*[n_points];
   std::complex<double> **wf_st_1=new std::complex<double>*[grid_size];
   std::complex<double> **wf_st_2=new std::complex<double>*[grid_size];
   pice_set *pice_data;

   for(int i=0;i!=grid_size;i++)
   {
      wf_st_1[i]=new std::complex<double>[n_times];
      wf_st_2[i]=new std::complex<double>[n_times];
   }
   for(int i=0;i!=n_theta;i++)
   {
      for(int j=0;j!=n_phi;j++)
      {
          thet[i*n_phi+j]=i*Pi/n_theta;
          phi[i*n_phi+j]=j*2*Pi/n_phi;
      }
   }
   for(int i=0;i!=n_phi;i++)
   {
   }
   for(int i=0;i!=n_points;i++)
   {
      dipole_1_x[i]=new std::complex<double>[grid_size];
      dipole_1_y[i]=new std::complex<double>[grid_size];
      dipole_1_z[i]=new std::complex<double>[grid_size];
      dipole_2_x[i]=new std::complex<double>[grid_size];
      dipole_2_y[i]=new std::complex<double>[grid_size];
      dipole_2_z[i]=new std::complex<double>[grid_size];
      mfpad[i]=new double[n_times];
   }
   pice_data= new pice_set(ionization_coupling_file);
   std::cout<<"GOT PICE DATA"<<std::endl;

   /*
   int x(0);
#pragma omp parallel for private(x)
   for( x=0;x<grid_size;x++)
   {
      for(int i=0;i<n_points;i++)
      {
         pice_data->fill_pice(&dipole_1_x[i][x],&dipole_1_y[i][x],&dipole_1_z[i][x],x,state_index_1,state_index_cat,thet[i],phi[i],kp,NULL);
         pice_data->fill_pice(&dipole_2_x[i][x],&dipole_2_y[i][x],&dipole_2_z[i][x],x,state_index_2,state_index_cat,thet[i],phi[i],kp,NULL);
      }
      std::cout<<"pice position "<<x<<" filled"<<std::endl;
   }
   std::cout<<"ALL PICE FILLED"<<std::endl;
*/
//EXTRACT THE NUCLEAR WAVE FUNCTION ON STATE 1
   input.open(state_1_wf_address.c_str());
   if(!input.is_open())
   {
      std::cout<<"FILE 1 ERROR"<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int t=0;t!=n_times;t++)
      {
         norm1[t]=0;
         for(int r=0;r!=grid_size;r++)
         {
            input>>temp;
            input>>temp;
            input>>rectemp;
            input>>imctemp;
            wf_st_1[r][t]=std::complex<double>(rectemp,imctemp);
            norm1[t]+=std::norm(wf_st_1[r][t]);
         }
      }
      input.close();
   }
   
   std::cout<<"STATE 1 WF DONE"<<std::endl;
//EXTRACT THE NUCLEAR WAVE FUNCTION ON STATE 2
   input.open(state_2_wf_address.c_str());
   if(!input.is_open())
   {
      std::cout<<"FILE 2 ERROR"<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int t=0;t!=n_times;t++)
      {
         norm2[t]=0;
         for(int r=0;r!=grid_size;r++)
         {
            input>>temp;
            input>>temp;
            input>>rectemp;
            input>>imctemp;
            wf_st_2[r][t]=std::complex<double>(rectemp,imctemp);
            norm2[t]+=std::norm(wf_st_2[r][t]);
         }
      }
      input.close();
   }
   std::cout<<"STATE 2 WF DONE"<<std::endl;

   std::cout<<"OUTPUT OF OVERLAP AS A FUNCTION OF TIME"<<std::endl;
   for(int t=0;t!=n_times;t++)
   {
      overlap=0;
      for(int r=0;r!=grid_size;r++)
      {
         overlap+=sqrt(std::norm(wf_st_1[r][t])*(std::norm(wf_st_2[r][t]))/(norm1[t]*norm2[t]));
      }
      std::cout<<0.04*t<<","<<overlap<<std::endl;
   }

   return 0;

   output.open(tdmfpad_address.c_str());
   output2.open(anisotropy_address.c_str());
   for(int t=0;t!=n_times;t++)
   {

      up=0;
      down=0;
      std::cout<<"WRITING OUTPUT TIME "<<t<<" / "<<n_times<<std::endl;
      for(int i=0;i!=n_points;i++)
      {
         mfpad[i][t]=0;
         for(int r=0;r!=grid_size;r++)
         {
            mfpad[i][t]+=std::norm(wf_st_1[r][t]*dipole_1_z[i][r])+std::norm(wf_st_2[r][t]*dipole_2_z[i][r])+2*std::real(wf_st_1[r][t]*dipole_1_z[i][r]*std::conj(wf_st_2[r][t]*dipole_2_z[i][r]));
         }
         if(i%n_phi==0 )
            output<<std::endl;
         output<<t*0.04<<"   "<<thet[i]<<"   "<<phi[i]<<"   "<<mfpad[i][t]<<std::endl;

         if(i<n_points/2)
         {
            up+=mfpad[i][t];
         }
         else
         {
            down+=mfpad[i][t];
         }
      }output<<std::endl;
      output2<<t*0.04<<"   "<<(up-down)/(up+down)<<std::endl;
   }
   output.close();
   output2.close();
//   std::cout<<"Electric field components : (X-comp) , (Y-comp) , (Z-comp)"<<std::endl<<x_comp<<" , "<<y_comp<<", "<<z_comp<<std::endl;


//        compute_mfpad(x_comp,y_comp,z_comp,n_points,mfpad,kp,dipole_address);//Compute the angularly resolved differential cross section for the given direction of the electric field

   return 0;
}

bool angularly_resolved_dipole_reader(double *prefactor,double *theta,double *phi, double * redipole_x,double * imdipole_x,double * redipole_y,double * imdipole_y,double * redipole_z,double * imdipole_z,int n_points,std::string dipole_address)
{
   using namespace std;
   string tmp_str;

   ifstream dipole_file;
   dipole_file.open(dipole_address.c_str());
   if(!dipole_file.is_open())
   {
      std::cout<<"CANNOT OPEN ANGULARLY RESOLVED IONIZATION DIPOLE FILE"<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_FAILURE);
   }
   for(int i=0;i!=n_points;i++)
   {
      dipole_file>>tmp_str;
      dipole_file>>theta[i];
      dipole_file>>phi[i];
      if(i==0)
         dipole_file>>*prefactor;
      else
         dipole_file>>tmp_str;
      dipole_file>>redipole_x[i];
      dipole_file>>imdipole_x[i];
      dipole_file>>redipole_y[i];
      dipole_file>>imdipole_y[i];
      dipole_file>>redipole_z[i];
      dipole_file>>imdipole_z[i];
   }
   dipole_file.close();
   return 1;
   
}
bool compute_mfpad(double x_comp,double y_comp,double z_comp,int n_points,double *mfpad,int kp,std::string dipole_address)
{
   int nk(50);
   int n_theta(20);
   int n_phi(40);
   double *RePICE_x=new double [n_points];
   double *ImPICE_x=new double [n_points];
   double *RePICE_y=new double [n_points];
   double *ImPICE_y=new double [n_points];
   double *RePICE_z=new double [n_points];
   double *ImPICE_z=new double [n_points];
   double *k=new double [n_points];
   double *thet=new double [n_points];
   double *phi=new double [n_points];

   using namespace std;
   ifstream input;
   input.open(dipole_address.c_str());

   if(!input.is_open())
   {      
      std::cout<<"CANNOT OPEN ANGULARLY RESOLVED IONIZATION DIPOLE FILE"<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int i=0;i!=nk*n_theta*n_phi;i++)
      {
         input>>k[i];
//         std::cout<<i<<"-"<<i%(n_theta*n_phi)<<"---"<<k[i]<<std::endl;
/*         if(i%(n_theta*n_phi)==0)
         {
            std::cout<<i/(n_theta*n_phi)<<"=>"<<k[i]<<std::endl;
         }*/
         input>>thet[i];
         input>>phi[i];
         input>>RePICE_x[i];
         input>>ImPICE_x[i];
         input>>RePICE_y[i];
         input>>ImPICE_y[i];
         input>>RePICE_z[i];
         input>>ImPICE_z[i];
      }
   input.close();
   for(int i=0;i!=n_theta*n_phi;i++)
   {
      mfpad[i]=real(std::complex<double>(x_comp*RePICE_x[kp*n_theta*n_phi+i]+y_comp*RePICE_y[kp*n_theta*n_phi+i]+z_comp*RePICE_z[kp*n_theta*n_phi+i],x_comp*ImPICE_x[kp*n_theta*n_phi+i]+y_comp*ImPICE_y[kp*n_theta*n_phi+i]+z_comp*ImPICE_z[kp*n_theta*n_phi+i])*std::complex<double>(x_comp*RePICE_x[kp*n_theta*n_phi+i]+y_comp*RePICE_y[kp*n_theta*n_phi+i]+z_comp*RePICE_z[kp*n_theta*n_phi+i],-(x_comp*ImPICE_x[kp*n_theta*n_phi+i]+y_comp*ImPICE_y[kp*n_theta*n_phi+i]+z_comp*ImPICE_z[kp*n_theta*n_phi+i])));
      std::cout<<k[kp*n_theta*n_phi+i]<<"   "<<thet[kp*n_theta*n_phi+i]<<"   "<<phi[kp*n_theta*n_phi+i]<<"   "<<mfpad[i]<<std::endl;
      if(i%n_phi==0 && i!=0)
         std::cout<<std::endl;
   }
   return 1;
}
}
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address)
{
   int n_theta(20);
   int n_phi(40);
   int temp_int;
   double * k=new double[n_points];
   double * theta=new double[n_points];
   double * phi=new double[n_points];
   double *RePICE_x=new double[n_points];
   double *ImPICE_x=new double[n_points];
   double *RePICE_y=new double[n_points];
   double *ImPICE_y=new double[n_points];
   double *RePICE_z=new double[n_points];
   double *ImPICE_z=new double[n_points];

   const double pi(acos(-1));
   const int nk(50);
   double total(0);
   double *cs=new double[nk];

   using namespace std;
   ifstream input;
   input.open("/data1/home/stephan/kxkykz.txt");
    double *distrib_cart=new double[3];

    for(int t=0;t!=2000;t++)
    {
       input>>temp_int;
       input>>distrib_cart[0];
       input>>distrib_cart[1];
       input>>distrib_cart[2];

       for(int p=0;p!=nk;p++)
       {
          theta[p*2000+t]=acos(distrib_cart[2]);
          phi[p*2000+t]=atan2(distrib_cart[1],distrib_cart[0]);
          if(phi[p*2000+t]<0)
             phi[p*2000+t]+=2*acos(-1);
       }
    }
    input.close();
    for(int p=0;p!=nk;p++)
    {
       for(int t=0;t!=2000;t++)
       {
          k[p*2000+t]=p*1.5/nk;
       }
    }

   input.open(dipole_address.c_str());

   if(!input.is_open())
   {      
      std::cout<<"CANNOT OPEN ANGULARLY RESOLVED IONIZATION DIPOLE FILE"<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int i=0;i!=nk*2000;i++)
      {
//         input>>k[i];
//         std::cout<<i<<"-"<<i%(n_theta*n_phi)<<"---"<<k[i]<<std::endl;
         if(i%(2000)==0)
         {
            std::cout<<i/(2000)<<"=>"<<k[i]<<std::endl;
         }
//         input>>theta[i];
//         input>>phi[i];
         input>>RePICE_x[i];
         input>>ImPICE_x[i];
         input>>RePICE_y[i];
         input>>ImPICE_y[i];
         input>>RePICE_z[i];
         input>>ImPICE_z[i];
      }
      for(int i=0;i!=nk;i++)
      {
         cs[i]=0;
         for(int j=0;j!=2000;j++)
         {
//             for(int l=0;l!=n_phi;l++)
             {
//                std::cout<<k[i*n_theta*n_phi+j*n_phi+l]<<","<<theta[i*n_theta*n_phi+j*n_phi+l]<<","<<phi[i*n_theta*n_phi+j*n_phi+l]<<","<<std::endl;
      //          cs[i]+=pow(k[i*n_theta*n_phi+j*n_phi+l],2)*sin(theta[i*n_theta*n_phi+j*n_phi+l])*(acos(-1)/n_theta)*(2*acos(-1)/n_phi)*real(std::complex<double>(x_comp*RePICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*RePICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*RePICE_z[i*n_theta*n_phi+j*n_phi+l],x_comp*ImPICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*ImPICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*ImPICE_z[i*n_theta*n_phi+j*n_phi+l])*std::complex<double>(x_comp*RePICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*RePICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*RePICE_z[i*n_theta*n_phi+j*n_phi+l],-(x_comp*ImPICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*ImPICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*ImPICE_z[i*n_theta*n_phi+j*n_phi+l])));
                cs[i]+=pow(k[i*2000+j],2)*(4*acos(-1)/2000)*real(std::complex<double>(x_comp*RePICE_x[i*2000+j]+y_comp*RePICE_y[i*2000+j]+z_comp*RePICE_z[i*2000+j],x_comp*ImPICE_x[i*2000+j]+y_comp*ImPICE_y[i*2000+j]+z_comp*ImPICE_z[i*2000+j])*std::complex<double>(x_comp*RePICE_x[i*2000+j]+y_comp*RePICE_y[i*2000+j]+z_comp*RePICE_z[i*2000+j],-(x_comp*ImPICE_x[i*2000+j]+y_comp*ImPICE_y[i*2000+j]+z_comp*ImPICE_z[i*2000+j])));
             }
         }
//         total+=cs[i]*(k[nk*n_theta*n_phi-1]-k[0])/nk;
         total+=cs[i]*(k[nk*2000-1]-k[0])/nk;
//         std::cout<<k[i*2000]<<", "<<cs[i]<<" => total = "<<total<<std::endl;
         std::cout<<k[i*2000]*k[i*2000]*27.211/2<<", "<<cs[i]<<std::endl;
 //        std::cout<<pow(k[i*n_theta*n_phi],2)*27.211/2<<", "<<cs[i]<<std::endl;
         //std::cout<<pow(k[i*n_theta*n_phi],2)*27.211/2<<", "<<cs[i]*Lx*Ly*Lz/pow((2*pi),3)<<std::endl;
      }
      input.close(); 
   }
   return 1;
}
double total_cs(int n_theta,int n_phi,double *theta,double *phi,double *mfpad)
{
   double sum(0);
   double Pi=acos(-1);
   for(int i=0;i!=n_theta;i++)
   {
      for(int j=0;j!=n_phi;j++)
      {
          sum+=sin(theta[i*n_phi+j])*(Pi/n_theta)*(2*Pi/n_phi)*mfpad[i*n_phi+j];
      }
   }
   return sum;
}
