#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>
#include "mkl.h"
#include "wf_overlap_headers.hpp"

void mo_cube_overlap_matrix(int n_sym,int *n_occs,double *mo_overlap,std::string cube_set_1,std::string cube_set_2);
void slater_det_overlap(int n_states,int n_occ,int ci_size_1,int ci_size_2,int n_elec,double **ci_vec_1,double **ci_vec_2, double *mo_overlap,double *det_overlap);
double wf_overlap_phase(std::string molpro_output_1,std::string molpro_output_2,std::string molpro_ovlp,int n_sym,int n_elec,int *n_states_sym,bool phase);

int main(int argc,char* argv[])
{
   int n_arg_min(11);
   //#######MOLPRO AND CUBE FILE NAMES AND PATHS#######
   std::string molpro_output_1("");
   std::string molpro_output_2("");
   std::string molpro_ovlp("");
   std::string molpro_cube_1("");
   std::string molpro_cube_2("");
   //#######MOLPRO COMPUTATION PARAMETERS AND VARIABLES#######
   int n_sym(0);
   int n_elec(0);
   int *n_states_sym;
   int n_states(0);
   //#######CHECK ARGC AND AGRV ARGUMENTS
   std::string temp;
   bool phase(0);
   bool n_sym_check(0);
   bool n_states_check(0);
   bool n_elec_check(0);
   bool molpro_out_1_check(0);
   bool molpro_out_2_check(0);
   bool molpro_cube_1_check(0);
   bool molpro_cube_2_check(0);
   bool molpro_ovlp_check(0);

   if(argc < n_arg_min)
   {
      std::cout<<"OVERLAP PROGRAM CALLING ERROR: NOT ENOUGH ARGUMENTS IN CALL"<<" "<<argc<<" / "<<n_arg_min<<std::endl;
      for(int i=0;i!=argc;i++)
      {
         std::cout<<"  "<<argv[i]<<std::endl;
      }
      exit(EXIT_FAILURE);
   }
   else
   {
      //std::cout<<"OVERLAP PROGRAM CALLED WITH "<<std::endl<<" "<<argc<<" ARGUMENTS "<<std::endl;

      for(int i=1;i!=argc;i++)
      {
         temp=argv[i];
         //std::cout<<"  "<<temp<<std::endl;
         if(temp == "n_sym")
         {
        //    std::cout<<"n_sym = "<<argv[i+1]<<std::endl;
            n_sym=atoi(argv[i+1]);
            i++;
            n_sym_check=1;
         }
         else if(temp == "n_elec")
         {
          //  std::cout<<"n_elec = "<<argv[i+1]<<std::endl;
            n_elec=atoi(argv[i+1]);
            i++;
            n_elec_check=1;
         }
         else if(temp == "n_states")
         {
            if(!n_sym_check)
            {
               std::cout<<"OVERLAP PROGRAM CALLING ERROR: N_SYM SHOULD BE GIVEN BEFORE N_STATES IN CALL"<<std::endl;
               exit(EXIT_FAILURE);
            }
            else
            {
               n_states_sym=new int[n_sym];
               for(int s=0;s!=n_sym;s++)
               {
            //       std::cout<<"n_states sym "<<s<<" = "<<argv[i+1+s]<<std::endl;
                   n_states_sym[s]=atoi(argv[i+1+s]);
               }
               i+=n_sym;
               n_states_check=1;
            }
         }
         else if(temp == "molpro_output_1")
         {
           // std::cout<<"molpro_output_1 = "<<argv[i+1]<<std::endl;
            molpro_output_1 = argv[i+1];
            molpro_out_1_check=1;
            i++;
         }
         else if(temp == "molpro_output_2")
         {
          //  std::cout<<"molpro_output_2 = "<<argv[i+1]<<std::endl;
            molpro_output_2 = argv[i+1];
            molpro_out_2_check=1;
            i++;
         }
         else if(temp == "molpro_ovlp")
         {
         //   std::cout<<"molpro_ovlp = "<<argv[i+1]<<std::endl;
            molpro_ovlp = argv[i+1];
            molpro_ovlp_check=1;
            i++;
         }
         else if(temp == "molpro_cube_1")
         {
            std::stringstream sscube1;
         //   std::cout<<"molpro_cube_1 = "<<argv[i+1]<<std::endl;
            sscube1.str("");
            sscube1<<argv[i+1]<<"_orbital_";
            molpro_cube_1=sscube1.str();
            molpro_cube_1_check=1;
            i++;
         }
         else if(temp == "molpro_cube_2")
         {
            std::stringstream sscube2;
        //    std::cout<<"molpro_cube_2 = "<<argv[i+1]<<std::endl;
            sscube2.str("");
            sscube2<<argv[i+1]<<"_orbital_";
            molpro_cube_2=sscube2.str(); 
            molpro_cube_2_check=1;
            i++;
         }
         else if(temp == "phase")
         {
         //   std::cout<<"Compute the phase between the electronic states. Will return either -1 or +1"<<std::endl;
            phase=1;
         }
         else
         {
            std::cout<<"OVERLAP PROGRAM CALLING ERROR: UNKNOWN PARAMETER : "<<temp<<std::endl;
            exit(EXIT_FAILURE);
         }
      }
   }

   wf_overlap_phase(molpro_output_1,molpro_output_2,molpro_ovlp,n_sym,n_elec,n_states_sym,phase);
   return 0;
}

double wf_overlap_phase(std::string molpro_output_1,std::string molpro_output_2,std::string molpro_ovlp,int n_sym,int n_elec,int *n_states_sym, bool phase)
{
   //########################################################################################
   //THIS ROUTINE COMPUTES THE OVERLAP BETWEEN ONE ELECTRONIC STATE COMPUTED IN TWO DIFFERENT MOLPRO FILES (E.G. : TWO DIFFERENT GEOMETRIES; COMPUTATION LEVELS,...)
   //WE COMPUTE THE OVERLAP <PSI(1)|PSI(2)>=∑(_K^CI)∑(_L^CI) C*(K).C(L)<D(1,K)|D(2,L)> WHERE D(X,K) IS A SLATER DETERMINANT. SINCE THE OVERLAP BEWTEEN SLATER DETERMINANTS IS THE DETERMINANT OF THE OVERLAPS BETWEEN MOLECULAR ORBITALS, WE COMPUTE THE OVERLAP BETWEEN MO'S BY USING CUBE FILES.
   //
   //                |<PHI(1;1,K)|PHI(1;2,L)>    <PHI(1;1,K)|PHI(2;2,L)>    ...........    <PHI(1;1,K)|PHI(N;2,L)>|
   //                |<PHI(2;1,K)|PHI(1;2,L)>    <PHI(2;1,K)|PHI(2;2,L)>    ...........    <PHI(2;1,K)|PHI(N;2,L)>|
   //                |           .                          .               .                                     |
   //                |           .                          .                 .                                   |
   //<D(1,K)|D(2,L)>=|           .                          .                   .                                 |
   //                |           .                          .                     .                               |
   //                |           .                          .                       .                             |
   //                |           .                          .                         .                           |
   //                |<PHI(N;1,K)|PHI(1;2,L)>    <PHI(N;1,K)|PHI(2;2,L)>    ...........    <PHI(N;1,K)|PHI(N;2,L)>|
   //
   //########################################################################################
   //std::cout<<"Entering overlap routine"<<std::endl;

   using namespace std;

   //#######MOLPRO COMPUTATION PARAMETERS AND VARIABLES#######
   int n_states;
   int *n_occs=new int [n_sym];
   int n_occ(0);
   int *n_closeds = new int[n_sym];
   int n_closed(0);
   int basis_size(0);
   int* num_of_ci_sym_1= new int[n_sym];
   int* num_of_ci_sym_2= new int[n_sym];
   int num_of_ci_1(0);
   int num_of_ci_2(0);
   double *ci_vector_1[2];
   double *ci_vector_2[2];
   double *mo_overlap;
   double *det_overlap;
   double *ES_overlap;

   if(n_sym == 1)
   {
      n_states = *n_states_sym;
   }
   else
   {
      n_states = 0;
      for(int s=0;s!=n_sym;s++)
      {
         n_states += n_states_sym[s];
      }
   }

   ES_overlap=new double[n_states*n_states];

   //std::cout<<"Array initialized"<<std::endl;
   //GET THE NUMBER OF CLOSED AND OCCUPIED MO'S
   size_query(n_occs,n_closeds,&basis_size,molpro_output_1,n_sym);
   for(int i=0;i!=n_sym;i++)
   {
      n_occ+=n_occs[i];
      n_closed+=n_closeds[i];
   }
   mo_overlap=new double[n_occ*n_occ];
   
   //std::cout<<"basis set size = "<<basis_size<<std::endl<<"n_closed = "<<n_closed<<std::endl<<"n_occ = "<<n_occ<<std::endl;
   //GET THE SIZE OF THE CI VECTORS AND INITIALIZE CI VECTORS
   num_of_ci_reader(n_states_sym,NULL,num_of_ci_sym_1,NULL,molpro_output_1,n_occs,n_sym);
   num_of_ci_reader(n_states_sym,NULL,num_of_ci_sym_2,NULL,molpro_output_2,n_occs,n_sym);

   //std::cout<<"Number of CI obtained"<<std::endl;

   for(int i=0;i!=n_sym;i++)
   {
      num_of_ci_1+=num_of_ci_sym_1[i];
      num_of_ci_2+=num_of_ci_sym_2[i];
      //std::cout<<"Geom1 : "<<num_of_ci_sym_1[i]<<" CI coeff. geom 2: "<<num_of_ci_sym_2[i]<<" CI coeff"<<std::endl;
   }
      ci_vector_1[0]=new double[n_elec*num_of_ci_1+n_states*num_of_ci_1];
      ci_vector_1[1]=new double[n_elec*num_of_ci_1];
      ci_vector_2[0]=new double[n_elec*num_of_ci_2+n_states*num_of_ci_2];
      ci_vector_2[1]=new double[n_elec*num_of_ci_2];
   for(int i=0;i!=n_elec*num_of_ci_1+n_states*num_of_ci_1;i++)
   {
      ci_vector_1[0][i]=0;
      ci_vector_2[0][i]=0;
   }
   det_overlap=new double[num_of_ci_1*num_of_ci_2];

   //std::cout<<"Initialized ci-vector size dependent arrays."<<std::endl;

   //GET THE CI COEFFICIENTS AND THE CONFIGURATIONS
   ci_vec_reader(n_states_sym,NULL,n_occs,n_closeds,n_elec,num_of_ci_sym_1,NULL,ci_vector_1,NULL,molpro_output_1,n_sym);
   ci_vec_reader(n_states_sym,NULL,n_occs,n_closeds,n_elec,num_of_ci_sym_2,NULL,ci_vector_2,NULL,molpro_output_2,n_sym);
   //std::cout<<"##################### "<<num_of_ci_1<<"CI Coefficients"<<std::endl;
   for(int m=0;m!=num_of_ci_1;m++)
   {
      for(int n=0;n!=n_states;n++)
      {
         std::cout<<setw(15)<<ci_vector_1[0][(n_elec+n_states)*m+n_elec+n]<<"   ";
         //std::cout<<ci_vector_1[0][(n_elec+n_states)*m+n]<<"   ";
      }std::cout<<"***"<<std::endl;
   }
   std::cout<<"#####################"<<num_of_ci_2<<"CI coefficients"<<std::endl;
   for(int m=0;m!=num_of_ci_2;m++)
   {
      for(int n=0;n!=n_states;n++)
      {
         std::cout<<setw(15)<<ci_vector_2[0][(n_elec+n_states)*m+n_elec+n]<<"   ";
         //std::cout<<ci_vector_2[0][(n_elec+n_states)*m+n]<<"   ";
      }std::cout<<"***"<<std::endl;
   }
   exit(EXIT_SUCCESS);
   //std::cout<<"#####################"<<std::endl;

   //std::cout<<"CI-vectors obtained"<<std::endl;
   //BUILD THE MO OVERLAP MATRIX FROM THE MO CUBES
   overlap_displaced_MO(mo_overlap,n_occs,&basis_size,molpro_output_1,molpro_output_2,molpro_ovlp,n_sym);
   //mo_cube_overlap_matrix(n_sym,n_occs,mo_overlap,molpro_cube_1,molpro_cube_2);
   //std::cout<<"MO overlap computed"<<std::endl;
   //std::cout<<"MO overlap matrix between the two geometries"<<std::endl;
/*   for(int i=0;i!=n_occ;i++)
   {
      for(int j=0;j!=n_occ;j++)
      {
         std::cout<<setw(15)<<setprecision(6)<<mo_overlap[n_occ*i+j];
         if(j%6 == 0 && j != 0)
            std::cout<<std::endl;

      }std::cout<<std::endl<<std::endl;
   }std::cout<<std::endl<<std::endl<<std::endl;*/
   
   //COMPUTE THE DETERMINANTS FOR EACH CONFIGURATION
   slater_det_overlap(n_states,n_occ,num_of_ci_1,num_of_ci_2,n_elec,ci_vector_1,ci_vector_2,mo_overlap,det_overlap);
   //std::cout<<"determinants overlap computed"<<std::endl;
   /*std::cout<<"Determinant overlap matrix between the two geometries"<<std::endl;
   for(int i=0;i!=num_of_ci_1;i++)
   {
      for(int j=0;j!=num_of_ci_2;j++)
      {
         std::cout<<setw(10)<<setprecision(6)<<det_overlap[num_of_ci_2*i+j];
      }std::cout<<std::endl;
   }std::cout<<std::endl<<std::endl;
   */
   //COMPUTE THE OVERLAP BETWEEN THE ELECTRONIC STATES
for(int m=0;m!=n_states;m++)
{
   for(int n=0;n!=n_states;n++)
   {
       ES_overlap[m*n_states+n]=0;
       for(int i=0;i!=num_of_ci_1;i++)
       {
          for(int j=0;j!=num_of_ci_2;j++)
          {
             //std::cout<<ci_vector_1[0][(n_elec+n_states)*i+n_elec+n]<<" ; "<<ci_vector_2[0][(n_elec+n_states)*j+n_elec+n]<<" ; "<<det_overlap[i*num_of_ci_2+j]<<std::endl;
             ES_overlap[m*n_states+n]+=ci_vector_1[0][(n_elec+n_states)*i+n_elec+m]*ci_vector_2[0][(n_elec+n_states)*j+n_elec+n]*det_overlap[i*num_of_ci_2+j];
          }
       }
   }
}

   if(!phase)
   {
      for(int i=0;i!=n_states;i++)
      {
         for(int j=0;j!=n_states;j++)
         {
            ES_overlap[j*n_states+i]=ES_overlap[i*n_states+j];
         }
      }
      for(int i=0;i!=n_states;i++)
      {
         for(int j=i;j!=n_states;j++)
         {
            //std::cout<<"OVERLAP BETWEEN THE STATES "<<j<<" AND "<<i<<" = "<<ES_overlap[i*n_states+j]<<std::endl; 
            std::cout<<setw(15)<<ES_overlap[i*n_states+j]<<std::endl; 
         }
      }
   }
   else
   {
      for(int i=0;i!=n_states;i++)
      {
         std::cout<<setw(15)<<ES_overlap[i*n_states+i]<<std::endl; 
      }
   }
}

void mo_cube_overlap_matrix(int n_sym,int *n_occs,double *mo_overlap,std::string cube_set_1,std::string cube_set_2)
{
   int num_of_atoms(0);
   int nx(0);
   int ny(0);
   int nz(0);
   double xmin(0);
   double ymin(0);
   double zmin(0);
   double xmax(0);
   double ymax(0);
   double zmax(0);
   double **mo_cube_set_1;
   double **mo_cube_set_2;
   int mo_index(0);
   int n_occ(0);

   using namespace std;

   for(int i=0;i!=n_sym;i++)
   {
      n_occ+=n_occs[i];
   }
   mo_cube_set_1=new double *[n_occ];
   mo_cube_set_2=new double *[n_occ];

   stringstream sizefile_name;
   string sizefile_name_str;
   sizefile_name.str("");
   sizefile_name<<cube_set_1.c_str()<<"1.1.cube";
   sizefile_name_str=sizefile_name.str();

   mo_cube_size_reader(&num_of_atoms,&xmin,&ymin,&zmin,&xmax,&ymax,&zmax,&nx,&ny,&nz,sizefile_name_str);
   std::cout<<"Cube size found. "<<std::endl<<"nx = "<<nx<<std::endl<<"ny = "<<ny<<std::endl<<"nz = "<<nz<<std::endl<<"xmin = "<<xmin<<std::endl<<"ymin = "<<ymin<<std::endl<<"zmin = "<<zmin<<std::endl<<"xmax = "<<xmax<<std::endl<<"ymax = "<<ymax<<std::endl<<"zmax = "<<zmax<<std::endl;
   for(int k=0;k!=n_occ;k++)
   {
      mo_cube_set_1[k]=new double[nx*ny*nz];
      mo_cube_set_2[k]=new double[nx*ny*nz];
   }
   mo_index=0;
   for(int i=0;i!=n_sym;i++)
   {
      for(int k=0;k!=n_occs[i];k++)
      {
         cube_reader(k,i,nx,ny,nz,cube_set_1,mo_cube_set_1[mo_index]);
//         std::cout<<mo_cube_set_1[mo_index][nx*ny*nz-1]<<"is the first element of geom 1 cube "<<mo_index<<std::endl;
         cube_reader(k,i,nx,ny,nz,cube_set_2,mo_cube_set_2[mo_index]);
//         std::cout<<mo_cube_set_2[mo_index][nx*ny*nz-1]<<"is the first element of geom 2 cube "<<mo_index<<std::endl;
         mo_index++;
      }
   }
   double norm1(0);
   double norm2(0);
   for(int i=0;i!=n_occ;i++)
   {
      cube_dot_product(mo_cube_set_1[i],mo_cube_set_1[i],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,1,&norm1);
//      std::cout<<"norm of MO "<<i<<" = "<<norm<<"(geom 1)"<<std::endl;
      for(int j=0;j!=n_occ;j++)
      {
         cube_dot_product(mo_cube_set_2[j],mo_cube_set_2[j],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,1,&norm2);
//         std::cout<<"norm of MO "<<i<<" = "<<norm<<"(geom 2)"<<std::endl;
//       cube_dot_product(mo_cube_set_2[i],mo_cube_set_2[j],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,1,&norm);
         cube_dot_product(mo_cube_set_1[i],mo_cube_set_2[j],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,1,&mo_overlap[n_occ*i+j]);
         mo_overlap[n_occ*i+j]/=(sqrt(norm1*norm2));
         std::cout<<mo_overlap[n_occ*i+j]<<"    ";
      }std::cout<<"***"<<std::endl;
   }

   delete [] mo_cube_set_1;
   delete [] mo_cube_set_2;
}
void slater_det_overlap(int n_states,int n_occ,int ci_size_1,int ci_size_2,int n_elec,double **ci_vec_1,double **ci_vec_2, double *mo_overlap,double *det_overlap)
{
   bool test(0);
   bool test2(0);


                
                for (int n=0; n!=ci_size_1; n++)//   over configurations of input 1
                {
                    for (int l=0; l!=ci_size_2; l++)//  over configuration of input 2
                    {
                       double *temp;
                       temp=new double[(n_elec)*(n_elec)];
                       det_overlap[ci_size_2*n+l]=0;
                        test2=0;
                     //========================VVVVVVVVVVVVV Overlap matrix for a given configuration VVVVVVVVVVVVVVVVV==================   
                        for (int m=0; m!=n_elec; m++)  //Over the electrons of input 1
                        {
                           for (int o=0; o!=n_elec; o++)//Over the electrons of input 2
                           {
                              temp[(n_elec)*m+o]=mo_overlap[n_occ*int(ci_vec_1[0][(n_elec+n_states)*n+m])+int(ci_vec_2[0][(n_elec+n_states)*l+o])]*kronecker_delta(ci_vec_1[1][n_elec*n+m], ci_vec_2[1][(n_elec)*l+o]);
                           }
                        }
                        //========================^^^^^^^^^^^^^^ Overlap matrix for a given configuration ^^^^^^^^^^^^^==================
                        delete [] temp;
                        
                        det_overlap[ci_size_2*n+l]=determinant(temp,(n_elec));
                    }
                }
}
/*
 * 
                        for (int m=0; m!=n_elec; m++)  //Over the electrons of input 1
                        {
                            for (int p=0; p!=n_occ; p++)//Over the MO of input 1
                            {
                                    if(int(ci_vec_1[0][(n_elec+n_states)*n+m])==p)
                                    {
                                        for (int o=0; o!=n_elec; o++)//Over the electrons of input 2
                                        {
                                            for (int q=0; q!=n_occ; q++)//Over the MO of input 2
                                            {
                                                if(int(ci_vec_2[0][(n_elec+n_states)*l+o])==q)
                                                {
                                                    temp[(n_elec)*m+o]=mo_overlap[n_occ*p+q]*kronecker_delta(ci_vec_1[1][n_elec*n+m], ci_vec_2[1][(n_elec)*l+o]);
                                                }
                                            }
                                        }
                                    }
                            }
                        }
 * */
