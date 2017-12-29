#include <cmath>


double hamilton_matrix::electric_field(double time_index,double* vector)
{
   double strength(0.015);
   double origin(10/0.02418884);
   double sigma(2/0.02418884);
   double energy(0.1);
   double CEP(0);

   double time(this->m_h*time_index);
   double amplitude(strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*(cos(energy*(time-origin)+CEP)-(time-origin)*sin(energy*(time-origin)+CEP)/(energy*sigma*sigma)));

   vector[0]=0;
   vector[1]=0;
   vector[2]=amplitude;

   return 0;
}
double hamilton_matrix::potential_vector(double time_index,double* vector)
{
   double strength(0.025);
   double origin(10/0.02418884);
   double sigma(2/0.02418884);
   double energy(0.1);
   double CEP(0);

   double time(this->m_h*time_index);
   double amplitude(strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*sin(energy*(time-origin)+CEP));
   vector[0]=0;
   vector[1]=0;
   vector[2]=amplitude;
   return 0;

}
