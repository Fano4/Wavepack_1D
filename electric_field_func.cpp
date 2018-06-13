#include <cmath>


double hamilton_matrix::electric_field(double time_index,double* vector)
{
   double strength(0.04);//0.022);
   double origin(3.4/0.02418884);//(4/0.02418884);//3.4/0.02418884);
   double sigma(0.677/0.02418884);//0.887/0.02418884);//0.677/0.02418884);
   double energy(0.05695);//0.08563);//0.05695);
   double CEP(acos(-1));

   double time(this->m_h*time_index);
   double amplitude(strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*(cos(energy*(time-origin)+CEP)-(time-origin)*sin(energy*(time-origin)+CEP)/(energy*sigma*sigma)));

   vector[0]=0;
   vector[1]=0;
   vector[2]=amplitude;

   return 0;
}
double hamilton_matrix::potential_vector(double time_index,double* vector)
{
   double strength(0.04);//0.022);
   double origin(3.4/0.02418884);//4/0.02418884);//3.4/0.02418884);
   double sigma(0.677/0.02418884);//0.887/0.02418884);//0.677/0.02418884);
   double energy(0.05695);//0.08563);//0.05695);
   double CEP(acos(-1));


   double time(this->m_h*time_index);
   double amplitude(strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*sin(energy*(time-origin)+CEP));
   vector[0]=0;
   vector[1]=0;
   vector[2]=amplitude;
   return 0;

}
