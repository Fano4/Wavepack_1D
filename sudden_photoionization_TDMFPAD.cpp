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
   string state_1_wf_address("/data1/home/stephan/wvpck_mfpad_probe/neut_test_pulse_2_wf_state_3.wvpck");
   string state_2_wf_address("/data1/home/stephan/wvpck_mfpad_probe/neut_test_pulse_2_wf_state_4.wvpck");
   string ionization_coupling_file("/data1/home/stephan/Wavepack_1D/wavepack_int_input/LiH_pice_data_newdyson.h5");//LiH_PICE_R_i_j.txt
   string tdmfpad_address("/data1/home/stephan/wvpck_mfpad_probe/tdmfpad_sudden_st_3_4.txt");
   string anisotropy_address("/data1/home/stephan/wvpck_mfpad_probe/anisotropy_sudden.txt");
   string tdspec_address("/data1/home/stephan/wvpck_mfpad_probe/tdpes_sudden_st_3_4.txt");

   ifstream input;
   ofstream output;
   ofstream output2;
   ofstream output3;

   //DIMENSION OF THE R GRID AND OF THE CONTINUUM
   int grid_size(512);
   int n_theta(20);
   int n_phi(40);
   int nk(50);
   double kmin(1e-4);
   double kmax(1.0);
   double kpp;
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
   std::complex<double> ***pesdipole_1_x=new std::complex<double>**[nk];
   std::complex<double> ***pesdipole_1_y=new std::complex<double>**[nk];
   std::complex<double> ***pesdipole_1_z=new std::complex<double>**[nk];
   std::complex<double> ***pesdipole_2_x=new std::complex<double>**[nk];
   std::complex<double> ***pesdipole_2_y=new std::complex<double>**[nk];
   std::complex<double> ***pesdipole_2_z=new std::complex<double>**[nk];
   double **tdpes=new double*[nk];
   double **mfpad=new double*[n_points];
   std::complex<double> **wf_st_1=new std::complex<double>*[grid_size];
   std::complex<double> **wf_st_2=new std::complex<double>*[grid_size];
   pice_set *pice_data;

   for(int i=0;i!=grid_size;i++)
   {
      wf_st_1[i]=new std::complex<double>[n_times];
      wf_st_2[i]=new std::complex<double>[n_times];
   }
   for(int i=0;i!=nk;i++)
   {
      tdpes[i]=new double[n_times];
      pesdipole_1_x[i]=new std::complex<double>*[n_points];
      pesdipole_1_y[i]=new std::complex<double>*[n_points];
      pesdipole_1_z[i]=new std::complex<double>*[n_points];
      pesdipole_2_x[i]=new std::complex<double>*[n_points];
      pesdipole_2_y[i]=new std::complex<double>*[n_points];
      pesdipole_2_z[i]=new std::complex<double>*[n_points];
      for(int j=0;j!=n_points;j++)
      {
         pesdipole_1_x[i][j]=new std::complex<double>[grid_size];
         pesdipole_1_y[i][j]=new std::complex<double>[grid_size];
         pesdipole_1_z[i][j]=new std::complex<double>[grid_size];
         pesdipole_2_x[i][j]=new std::complex<double>[grid_size];
         pesdipole_2_y[i][j]=new std::complex<double>[grid_size];
         pesdipole_2_z[i][j]=new std::complex<double>[grid_size];
      }
   }
   for(int i=0;i!=n_theta;i++)
   {
      for(int j=0;j!=n_phi;j++)
      {
          thet[i*n_phi+j]=i*Pi/n_theta;
          phi[i*n_phi+j]=j*2*Pi/n_phi;
      }
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

   
   int x(0);
#pragma omp parallel for private(x,kpp)
   for( x=0;x<grid_size;x++)
   {
      std::cout<<"pice position "<<x<<"...";
      for(int i=0;i<n_points;i++)
      {
//         std::cout<<" angle "<<i<<"...";
         pice_data->fill_pice(&dipole_1_x[i][x],&dipole_1_y[i][x],&dipole_1_z[i][x],x,state_index_1,state_index_cat,thet[i],phi[i],kp,NULL);
         pice_data->fill_pice(&dipole_2_x[i][x],&dipole_2_y[i][x],&dipole_2_z[i][x],x,state_index_2,state_index_cat,thet[i],phi[i],kp,NULL);
         for(int j=0;j<nk;j++)
         {
            kpp=kmin+j*(kmax-kmin)/nk;
            pice_data->fill_pice(&pesdipole_1_x[j][i][x],&pesdipole_1_y[j][i][x],&pesdipole_1_z[j][i][x],x,state_index_1,state_index_cat,thet[i],phi[i],kpp,NULL);
            pice_data->fill_pice(&pesdipole_2_x[j][i][x],&pesdipole_2_y[j][i][x],&pesdipole_2_z[j][i][x],x,state_index_2,state_index_cat,thet[i],phi[i],kpp,NULL);
//          std::cout<<"x = "<<x<<"; i = "<<i<<"; j = "<<j<<std::endl
//             <<"k = "<<kpp<<"; "<<pesdipole_1_z[j][i][x]<<" ; "<<pesdipole_2_z[j][i][x]<<std::endl;;
         }
      }
      std::cout<<"Done !"<<std::endl;
   }
   std::cout<<"ALL PICE FILLED"<<std::endl;

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

   output.open(tdmfpad_address.c_str());
   output2.open(anisotropy_address.c_str());
   output3.open(tdspec_address.c_str());
   for(int t=0;t!=n_times;t++)
   {

      up=0;
      down=0;
      std::cout<<"WRITING OUTPUT TIME "<<t<<" / "<<n_times<<std::endl;
      for(int i=0;i!=nk;i++)
      {
         kpp=kmin+i*(kmax-kmin)/nk;
         tdpes[i][t]=0;
         for(int r=0;r!=grid_size;r++)
         {
            for(int j=0;j!=n_points;j++)
            {
               std::cout<<"probe! i = "<<i<<" ; r = "<<r<<" ; j = "<<j
                  <<std::endl<<std::norm(pesdipole_1_z[i][j][r])<<";"<<std::norm(pesdipole_2_z[i][j][r])<<";"<<std::norm(wf_st_1[r][t])<<","<<std::norm(wf_st_2[r][t])<<std::endl;
               tdpes[i][t]+=(std::norm(wf_st_1[r][t]*pesdipole_1_z[i][j][r])+std::norm(wf_st_2[r][t]*pesdipole_2_z[i][j][r])+2*std::real(wf_st_1[r][t]*pesdipole_1_z[i][j][r]*std::conj(wf_st_2[r][t]*pesdipole_2_z[i][j][r])))*pow(kpp,2)*sin(thet[j])*2*pow(acos(-1),2)/(n_theta*n_phi);
            }
         }
         output3<<t*0.04<<"   "<<kpp<<"   "<<tdpes[i][t]<<std::endl;
      }output3<<std::endl;
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
   output3.close();


   return 0;
}
