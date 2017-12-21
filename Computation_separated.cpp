#include "Computation.hpp"

bool Runge_kutta(double **Rey0, double **Imy0, double h, double time, int* n_states, double **pot, double **dipole_x, double **dipole_y, double **dipole_z,double **RePICE_x,double **ImPICE_x,double **RePICE_y,double **ImPICE_y,double **RePICE_z,double **ImPICE_z,double *k,double *theta, double *phi,double *pulse1_param,double *pulse2_param)
{

   int n_states_neut(n_states[0]);
   int n_states_cat(n_states[1]);
   int n_states_cont(n_states[2]);

    double **k1=new double * [3];
    k1[0]=new double[n_states_neut];
    k1[1]=new double[n_states_cat];
    k1[2]=new double[n_states_cont];
    double **k2=new double * [3];
    k2[0]=new double[n_states_neut];
    k2[1]=new double[n_states_cat];
    k2[2]=new double[n_states_cont];
    double **k3=new double * [3];
    k3[0]=new double[n_states_neut];
    k3[1]=new double[n_states_cat];
    k3[2]=new double[n_states_cont];
    double **k4=new double * [3];
    k4[0]=new double[n_states_neut];
    k4[1]=new double[n_states_cat];
    k4[2]=new double[n_states_cont];


    double **k1p=new double * [3];
    k1p[0]=new double[n_states_neut];
    k1p[1]=new double[n_states_cat];
    k1p[2]=new double[n_states_cont];
    double **k2p=new double * [3];
    k2p[0]=new double[n_states_neut];
    k2p[1]=new double[n_states_cat];
    k2p[2]=new double[n_states_cont];
    double **k3p=new double * [3];
    k3p[0]=new double[n_states_neut];
    k3p[1]=new double[n_states_cat];
    k3p[2]=new double[n_states_cont];
    double **k4p=new double * [3];
    k4p[0]=new double[n_states_neut];
    k4p[1]=new double[n_states_cat];
    k4p[2]=new double[n_states_cont];

    double **temp1=new double * [3];
    temp1[0]=new double[n_states_neut];
    temp1[1]=new double[n_states_cat];
    temp1[2]=new double[n_states_cont];
    double **temp2=new double * [3];
    temp2[0]=new double[n_states_neut];
    temp2[1]=new double[n_states_cat];
    temp2[2]=new double[n_states_cont];

    double Rey1;
    double Imy1;

    const double a21(.5),a32(.5),a43(1.),ddc1(1./6.),ddc2(1./3.),ddc3(1./3.),ddc4(1./6.);

    t_deriv(Rey0,Imy0, k1p, k1,time,n_states,pot,dipole_x,dipole_y,dipole_z,RePICE_x,ImPICE_x,RePICE_y,ImPICE_y,RePICE_z,ImPICE_z,k,theta,phi,pulse1_param,pulse2_param);

    for (int n=0; n!=n_states_neut; n++)
    {
                temp1[0][n]=Imy0[0][n]+a21*h*k1[0][n];
                temp2[0][n]=Rey0[0][n]+a21*h*k1p[0][n];
    }
    for (int n=0; n!=n_states_cat; n++)
    {
                temp1[1][n]=Imy0[1][n]+a21*h*k1[1][n];
                temp2[1][n]=Rey0[1][n]+a21*h*k1p[1][n];
    }
    for (int n=0; n!=n_states_cont; n++)
    {
                temp1[2][n]=Imy0[2][n]+a21*h*k1[2][n];
                temp2[2][n]=Rey0[2][n]+a21*h*k1p[2][n];
    }

     t_deriv(temp2,temp1, k2p, k2,time+0.5*h,n_states,pot,dipole_x,dipole_y,dipole_z,RePICE_x,ImPICE_x,RePICE_y,ImPICE_y,RePICE_z,ImPICE_z,k,theta,phi,pulse1_param,pulse2_param);

     for (int n=0; n!=n_states_neut; n++)
     {
                 temp1[0][n]=Imy0[0][n]+a32*h*k2[0][n];
                 temp2[0][n]=Rey0[0][n]+a32*h*k2p[0][n];
     }
     for (int n=0; n!=n_states_cat; n++)
     {
                 temp1[1][n]=Imy0[1][n]+a32*h*k2[1][n];
                 temp2[1][n]=Rey0[1][n]+a32*h*k2p[1][n];
     }
     for (int n=0; n!=n_states_cont; n++)
     {
                 temp1[2][n]=Imy0[2][n]+a32*h*k2[2][n];
                 temp2[2][n]=Rey0[2][n]+a32*h*k2p[2][n];
     }

     t_deriv(temp2,temp1, k3p, k3,time+0.5*h,n_states,pot,dipole_x,dipole_y,dipole_z,RePICE_x,ImPICE_x,RePICE_y,ImPICE_y,RePICE_z,ImPICE_z,k,theta,phi,pulse1_param,pulse2_param);

     for (int n=0; n!=n_states_neut; n++)
     {
                 temp1[0][n]=Imy0[0][n]+a43*h*k3[0][n];
                 temp2[0][n]=Rey0[0][n]+a43*h*k3p[0][n];
     }
     for (int n=0; n!=n_states_cat; n++)
     {
                 temp1[1][n]=Imy0[1][n]+a43*h*k3[1][n];
                 temp2[1][n]=Rey0[1][n]+a43*h*k3p[1][n];
     }
     for (int n=0; n!=n_states_cont; n++)
     {
                 temp1[2][n]=Imy0[2][n]+a43*h*k3[2][n];
                 temp2[2][n]=Rey0[2][n]+a43*h*k3p[2][n];
     }

     t_deriv(temp2,temp1, k4p, k4,time+h,n_states,pot,dipole_x,dipole_y,dipole_z,RePICE_x,ImPICE_x,RePICE_y,ImPICE_y,RePICE_z,ImPICE_z,k,theta,phi,pulse1_param,pulse2_param);

     for (int n=0; n!=n_states_neut; n++)
     {
                     Imy1 = Imy0[0][n] + h*(ddc1*k1[0][n]+ddc2*k2[0][n]+ddc3*k3[0][n]+ddc4*k4[0][n]);
                     Rey1 = Rey0[0][n] + h*(ddc1*k1p[0][n]+ddc2*k2p[0][n]+ddc3*k3p[0][n]+ddc4*k4p[0][n]);

                     Imy0[0][n] = Imy1;
                     Rey0[0][n] = Rey1;
     }
     for (int n=0; n!=n_states_cat; n++)
     {
                     Imy1 = Imy0[1][n] + h*(ddc1*k1[1][n]+ddc2*k2[1][n]+ddc3*k3[1][n]+ddc4*k4[1][n]);
                     Rey1 = Rey0[1][n] + h*(ddc1*k1p[1][n]+ddc2*k2p[1][n]+ddc3*k3p[1][n]+ddc4*k4p[1][n]);

                     Imy0[1][n] = Imy1;
                     Rey0[1][n] = Rey1;
     }
     for (int n=0; n!=n_states_cont; n++)
     {
                     Imy1 = Imy0[2][n] + h*(ddc1*k1[2][n]+ddc2*k2[2][n]+ddc3*k3[2][n]+ddc4*k4[2][n]);
                     Rey1 = Rey0[2][n] + h*(ddc1*k1p[2][n]+ddc2*k2p[2][n]+ddc3*k3p[2][n]+ddc4*k4p[2][n]);

                     Imy0[2][n] = Imy1;
                     Rey0[2][n] = Rey1;
     }

    return 0;
}

bool t_deriv(double** Repsi, double** Impsi, double** Redpsi, double** Imdpsi, double time, int *n_states,double **pot, double **dipole_x, double **dipole_y, double **dipole_z,double **RePICE_x,double **ImPICE_x,double **RePICE_y,double **ImPICE_y,double **RePICE_z,double **ImPICE_z,double *k,double *theta, double *phi,double *pulse1_param,double *pulse2_param)
{
    double elec_field[3];
    double vector[3];
    int n_states_neut(n_states[0]);
    int n_states_cat(n_states[1]);
    int n_states_cont(n_states[2]);

    double **ImHpsi=new double * [3];
    ImHpsi[0]=new double [n_states_neut];
    ImHpsi[1]=new double [n_states_cat];
    ImHpsi[2]=new double [n_states_cont];
    double **ReHpsi=new double * [3];
    ReHpsi[0]=new double [n_states_neut];
    ReHpsi[1]=new double [n_states_cat];
    ReHpsi[2]=new double [n_states_cont];

    electric_field(elec_field, time,pulse1_param,pulse2_param);

    //QHQ
    for (int m=0;m!=n_states_neut;m++)
    {
        ImHpsi[0][m]=potential(m,0,pot)*Impsi[0][m];
        ReHpsi[0][m]=potential(m,0,pot)*Repsi[0][m];

        for (int n=0;n!=n_states_neut;n++)
        {
            ImHpsi[0][m] -= (elec_field[0]*dipole_moment(0,m,n,0,n_states,dipole_x,dipole_y,dipole_z) + elec_field[1]*dipole_moment(1,m,n,0,n_states,dipole_x,dipole_y,dipole_z) + elec_field[2]*dipole_moment(2,m,n,0,n_states,dipole_x,dipole_y,dipole_z))*Impsi[0][n];
            ReHpsi[0][m] -= (elec_field[0]*dipole_moment(0,m,n,0,n_states,dipole_x,dipole_y,dipole_z) + elec_field[1]*dipole_moment(1,m,n,0,n_states,dipole_x,dipole_y,dipole_z) + elec_field[2]*dipole_moment(2,m,n,0,n_states,dipole_x,dipole_y,dipole_z))*Repsi[0][n];
        }
    }



    //PHP Cation states
    for (int m=0;m!=n_states_cat;m++)
    {
        ImHpsi[1][m]=potential(m,1,pot)*Impsi[1][m];
        ReHpsi[1][m]=potential(m,1,pot)*Repsi[1][m];

        for (int n=0;n!=n_states_cat;n++)
        {
            ImHpsi[1][m] -= (elec_field[0]*dipole_moment(0,m,n,1,n_states,dipole_x,dipole_y,dipole_z) + elec_field[1]*dipole_moment(1,m,n,1,n_states,dipole_x,dipole_y,dipole_z) + elec_field[2]*dipole_moment(2,m,n,1,n_states,dipole_x,dipole_y,dipole_z))*Impsi[1][n];
            ReHpsi[1][m] -= (elec_field[0]*dipole_moment(0,m,n,1,n_states,dipole_x,dipole_y,dipole_z) + elec_field[1]*dipole_moment(1,m,n,1,n_states,dipole_x,dipole_y,dipole_z) + elec_field[2]*dipole_moment(2,m,n,1,n_states,dipole_x,dipole_y,dipole_z))*Repsi[1][n];
        }
    }


    //PHP Continuum states
    //potential_vector(vector,time,pulse1_param,pulse2_param);
    for (int m=0;m!=n_states_cont;m++)
    {
        ImHpsi[2][m]+=(pow(k[m]*sin(theta[m])*cos(phi[m])/*-vector[0]*/,2)+pow(k[m]*sin(theta[m])*sin(phi[m])/*-vector[1]*/,2)+pow(k[m]*cos(theta[m])/*-vector[2]*/,2))*Impsi[2][m];
        ReHpsi[2][m]+=(pow(k[m]*sin(theta[m])*cos(phi[m])/*-vector[0]*/,2)+pow(k[m]*sin(theta[m])*sin(phi[m])/*-vector[1]*/,2)+pow(k[m]*cos(theta[m])/*-vector[2]*/,2))*Repsi[2][m];
    }
    
    //PHQ Neutral to cation
    for(int m=0;m!=n_states_neut;m++)
    {
       for(int n=0;n!=n_states_cat;n++)
       {
          for(int o=0;o!=n_states_cont;o++)
          {
             ReHpsi[1][n]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Repsi[0][m];
             ImHpsi[1][n]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Impsi[0][m];
             ReHpsi[1][n]+=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Impsi[0][m];
             ImHpsi[1][n]-=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Repsi[0][m];
             ReHpsi[2][o]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Repsi[0][m];
             ImHpsi[2][o]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Impsi[0][m];
             ReHpsi[2][o]+=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Impsi[0][m];
             ImHpsi[2][o]-=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Repsi[0][m];
          }
       }
    }

    //QHP Cation to Neutral
    for(int m=0;m!=n_states_neut;m++)
    {
       for(int n=0;n!=n_states_cat;n++)
       {
          for(int o=0;o!=n_states_cont;o++)
          {
             ReHpsi[0][m]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Repsi[1][n]*Repsi[2][o];
             ReHpsi[0][m]-=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Repsi[1][n]*Impsi[2][o];
             ReHpsi[0][m]+=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Impsi[1][n]*Impsi[2][o];
             ReHpsi[0][m]-=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Impsi[1][n]*Repsi[2][o];
             ImHpsi[0][m]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Impsi[1][n]*Repsi[2][o];
             ImHpsi[0][m]-=(elec_field[0]*RePICE_x[m*n_states_cat+n][o] + elec_field[1]*RePICE_y[m*n_states_cat+n][o] + elec_field[2]*RePICE_z[m*n_states_cat+n][o])*Repsi[1][n]*Impsi[2][o];
             ImHpsi[0][m]+=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Repsi[1][n]*Repsi[2][o];
             ImHpsi[0][m]-=(elec_field[0]*ImPICE_x[m*n_states_cat+n][o] + elec_field[1]*ImPICE_y[m*n_states_cat+n][o] + elec_field[2]*ImPICE_z[m*n_states_cat+n][o])*Impsi[1][n]*Impsi[2][o];
          }
       }
    }
    for (int m=0; m!=n_states_neut;m++)
    {
        Imdpsi[0][m]=-ReHpsi[0][m];
        Redpsi[0][m]=ImHpsi[0][m];
    }
    for (int m=0; m!=n_states_cat;m++)
    {
        Imdpsi[1][m]=-ReHpsi[1][m];
        Redpsi[1][m]=ImHpsi[1][m];
    }
    for (int m=0; m!=n_states_cont;m++)
    {
        Imdpsi[2][m]=-ReHpsi[2][m];
        Redpsi[2][m]=ImHpsi[2][m];
    }

    delete[] ReHpsi[0];
    delete[] ReHpsi[1];
    delete[] ReHpsi[2];
    delete[] ReHpsi;
    delete[] ImHpsi[0];
    delete[] ImHpsi[1];
    delete[] ImHpsi[2];
    delete[] ImHpsi;
    return 0;
}

void electric_field(double *vector,double time,double *pulse1_param,double *pulse2_param)
{
   double h(0.00000001);
   double vector_im2[3];
   double vector_im1[3];
   double vector_i[3];
   double vector_ip1[3];
   double vector_ip2[3];

   potential_vector(vector_im2,time-2*h,pulse1_param,pulse2_param);
   potential_vector(vector_im1,time-h,pulse1_param,pulse2_param);
   potential_vector(vector_ip1,time-h,pulse1_param,pulse2_param);
   potential_vector(vector_ip2,time+2*h,pulse1_param,pulse2_param);


    vector[0]=-(1/(12*h))*(8*vector_ip1[0]-vector_ip2[0]-8*vector_im2[0]+vector_im2[0]);//composante X


    vector[1]=/*-(1/(12*h))*(8*vector_ip1[1]-vector_ip2[1]-8*vector_im2[1]+vector_im2[1]);*/sin(20*acos(-1)/180)*pulse1_param[0]*exp(-pow((time-pulse1_param[2]),2)/(2*pow(pulse1_param[3],2)))*(cos(pulse1_param[1]*(time-pulse1_param[2])+pulse1_param[4])-(time-pulse1_param[2])*sin(pulse1_param[1]*(time-pulse1_param[2])+pulse1_param[4])/(pulse1_param[1]*pulse1_param[3]*pulse1_param[3]));//composante Y


    vector[2]=/*-(1/(12*h))*(8*vector_ip1[2]-vector_ip2[2]-8*vector_im2[2]+vector_im2[2]);*/cos(20*acos(-1)/180)*pulse1_param[0]*exp(-pow((time-pulse1_param[2]),2)/(2*pow(pulse1_param[3],2)))*(cos(pulse1_param[1]*(time-pulse1_param[2])+pulse1_param[4])-(time-pulse1_param[2])*sin(pulse1_param[1]*(time-pulse1_param[2])+pulse1_param[4])/(pulse1_param[1]*pulse1_param[3]*pulse1_param[3]));//composante Z

}

double potential(int state,int species,double **pot)
{
    //double pot[8]={0,0.31891634,0.32940871,0.32942163,0.35437268,0.35531179,0.38837289,0.38888324};
    return pot[species][state];
}

double dipole_moment(int component_0_x_1_y_2_z,int matrix_row,int matrix_column,int species,int *n_states,double **dipole_x,double **dipole_y,double **dipole_z)
{
    int temp;
    int matrix_size(n_states[species]);
    //(matrix_size * matrix_row) + matrix_column – ((matrix_row * (matrix_row+1)) / 2)

    if(matrix_row>matrix_column)
    {
        temp=matrix_row;
        matrix_row=matrix_column;
        matrix_column=temp;
    }
                        switch (component_0_x_1_y_2_z)
                         {
                             case 0:
                             return dipole_x[species][(matrix_size * matrix_row) + matrix_column - int(double(matrix_row * (matrix_row+1)) / 2)];
                             break;
                             case 1:
                             return dipole_y[species][(matrix_size * matrix_row) + matrix_column - int(double(matrix_row * (matrix_row+1)) / 2)];
                             break;
                             case 2:
                             return dipole_z[species][(matrix_size * matrix_row) + matrix_column - int(double(matrix_row * (matrix_row+1)) / 2)];
                             break;
                         }
    return 0;
}

double vector_prod(double vector1[],double vector2[],int gsize)
{
    double sum(0.);

        for (int j=0; j!=gsize; j++)
        {
            sum+=vector1[j]*vector2[j];
        }

    return sum;
}

void mean_dip_mom(double **Repsi,double **Impsi, double *tDMX,double *tDMY,double *tDMZ,int *n_states, double **dipole_x, double **dipole_y, double **dipole_z,double **RePICE_x,double ** ImPICE_x,double **RePICE_y,double ** ImPICE_y,double **RePICE_z,double ** ImPICE_z)
{
    *tDMX=0;
    *tDMY=0;
    *tDMZ=0;
    int n_states_neut(n_states[0]);
    int n_states_cat(n_states[1]);
    int n_states_cont(n_states[2]);

    for (int i=0; i!=n_states_neut; i++)
    {
        *tDMX+=dipole_moment(0,i,i,0,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[0][i]*Repsi[0][i]+Impsi[0][i]*Impsi[0][i]);
        *tDMY+=dipole_moment(1,i,i,0,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[0][i]*Repsi[0][i]+Impsi[0][i]*Impsi[0][i]);
        *tDMZ+=dipole_moment(2,i,i,0,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[0][i]*Repsi[0][i]+Impsi[0][i]*Impsi[0][i]);
        for(int j=i+1;j!=n_states_neut;j++)
        {
            *tDMX+=2*dipole_moment(0,i,j,0,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[0][i]*Repsi[0][j]+Impsi[0][i]*Impsi[0][j]);
            *tDMY+=2*dipole_moment(1,i,j,0,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[0][i]*Repsi[0][j]+Impsi[0][i]*Impsi[0][j]);
            *tDMZ+=2*dipole_moment(2,i,j,0,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[0][i]*Repsi[0][j]+Impsi[0][i]*Impsi[0][j]);
        }
    }
    for (int i=0; i!=n_states_cat; i++)
    {
        *tDMX+=dipole_moment(0,i,i,1,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[1][i]*Repsi[1][i]+Impsi[1][i]*Impsi[1][i]);
        *tDMY+=dipole_moment(1,i,i,1,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[1][i]*Repsi[1][i]+Impsi[1][i]*Impsi[1][i]);
        *tDMZ+=dipole_moment(2,i,i,1,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[1][i]*Repsi[1][i]+Impsi[1][i]*Impsi[1][i]);
        for(int j=i+1;j!=n_states_cat;j++)
        {
            *tDMX+=2*dipole_moment(0,i,j,1,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[1][i]*Repsi[1][j]+Impsi[1][i]*Impsi[1][j]);
            *tDMY+=2*dipole_moment(1,i,j,1,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[1][i]*Repsi[1][j]+Impsi[1][i]*Impsi[1][j]);
            *tDMZ+=2*dipole_moment(2,i,j,1,n_states,dipole_x,dipole_y,dipole_z)*(Repsi[1][i]*Repsi[1][j]+Impsi[1][i]*Impsi[1][j]);
        }
    }
    for(int i=0;i!=n_states_neut;i++)
    {
       for(int j=0;j!=n_states_cat;j++)
       {
          for(int k=0;k!=n_states_cont;k++)
          {
             *tDMX+=2*(RePICE_x[i*n_states_cat+j][k]*Repsi[0][i]*Repsi[1][j]*Repsi[2][k]-RePICE_x[i*n_states_cat+j][k]*Repsi[0][i]*Impsi[1][j]*Impsi[2][k]+RePICE_x[i*n_states_cat+j][k]*Impsi[0][i]*Repsi[1][j]*Impsi[2][k]+RePICE_x[i*n_states_cat+j][k]*Impsi[0][i]*Impsi[1][j]*Repsi[2][k]+ImPICE_x[i*n_states_cat+j][k]*Repsi[0][i]*Repsi[1][j]*Impsi[2][k]+ImPICE_x[i*n_states_cat+j][k]*Repsi[0][i]*Impsi[1][j]*Repsi[2][k]);

             *tDMY+=2*(RePICE_y[i*n_states_cat+j][k]*Repsi[0][i]*Repsi[1][j]*Repsi[2][k]-RePICE_y[i*n_states_cat+j][k]*Repsi[0][i]*Impsi[1][j]*Impsi[2][k]+RePICE_y[i*n_states_cat+j][k]*Impsi[0][i]*Repsi[1][j]*Impsi[2][k]+RePICE_y[i*n_states_cat+j][k]*Impsi[0][i]*Impsi[1][j]*Repsi[2][k]+ImPICE_y[i*n_states_cat+j][k]*Repsi[0][i]*Repsi[1][j]*Impsi[2][k]+ImPICE_y[i*n_states_cat+j][k]*Repsi[0][i]*Impsi[1][j]*Repsi[2][k]);

             *tDMZ+=2*(RePICE_z[i*n_states_cat+j][k]*Repsi[0][i]*Repsi[1][j]*Repsi[2][k]-RePICE_z[i*n_states_cat+j][k]*Repsi[0][i]*Impsi[1][j]*Impsi[2][k]+RePICE_z[i*n_states_cat+j][k]*Impsi[0][i]*Repsi[1][j]*Impsi[2][k]+RePICE_z[i*n_states_cat+j][k]*Impsi[0][i]*Impsi[1][j]*Repsi[2][k]+ImPICE_z[i*n_states_cat+j][k]*Repsi[0][i]*Repsi[1][j]*Impsi[2][k]+ImPICE_z[i*n_states_cat+j][k]*Repsi[0][i]*Impsi[1][j]*Repsi[2][k]);
          }
       }
    }

}
void potential_vector(double * vector,double time,double *pulse1_param,double *pulse2_param)
{
   //PULSE PARAM VECTOR CONTAINS THE PARAMETERS OF THE GAUSSIAN PULSE
   //PULSE_PARAM[0] -> INTENSITY
   //PULSE_PARAM[1] -> ENERGY
   //PULSE_PARAM[2] -> CENTER TIME
   //PULSE_PARAM[3] -> STANDARD DEVIATION
   //PULSE_PARAM[4] -> CEP

   const double Pi=acos(-1);
   int phase(int((time-pulse1_param[2])*pulse1_param[1]/(2*Pi)));

   vector[0]=0;
   vector[1]=0;
   vector[2]=pulse1_param[0]*sin(pulse1_param[1]*(time-pulse1_param[2]) +pulse1_param[4]-2*Pi*phase)*exp(-pow(time-pulse1_param[2],2)/(2*pulse1_param[3]*pulse1_param[3]));
}
