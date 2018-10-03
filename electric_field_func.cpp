#include <cmath>


double hamilton_matrix::electric_field(double time_index,double* vector)
{
   double strength(0.04);//0.022);
   double origin(3.4/0.02418884);//(4/0.02418884);//3.4/0.02418884);
   double sigma(0.677/0.02418884);//0.887/0.02418884);//0.677/0.02418884);
   double energy(0.05695);//0.08563);//0.05695);
   double CEP(0-2*acos(-1)*int((energy*origin)/(2*acos(-1))));//acos(-1));
//   double CEP(acos(-1)-2*acos(-1)*int((energy*origin)/(2*acos(-1))));//acos(-1));

   double pump_probe_delay(8.0/0.02418884);

   double strength_probe(0.000);//0.022);
   double origin_probe(origin+pump_probe_delay);//(4/0.02418884);//3.4/0.02418884);
   double sigma_probe(0.677/0.02418884);//0.887/0.02418884);//0.677/0.02418884);
   double energy_probe(0.05695);//0.08563);//0.05695);
   double CEP_probe(0-2*acos(-1)*int((energy_probe*origin_probe)/(2*acos(-1))));//acos(-1));

   double time(this->m_h*time_index);
   double amplitude(strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*(cos(energy*(time-origin)+CEP)-(time-origin)*sin(energy*(time-origin)+CEP)/(energy*sigma*sigma))
          +strength_probe*exp(-pow((time-origin_probe),2)/(2*sigma_probe*sigma_probe))*(cos(energy_probe*(time-origin_probe)+CEP_probe)-(time-origin_probe)*sin(energy_probe*(time-origin_probe)+CEP_probe)/(energy_probe*sigma_probe*sigma_probe))
         );

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
   double energy(0.05695);//0.12656)/0.08563);//0.05695);
   double CEP(0-2*acos(-1)*int((energy*origin)/(2*acos(-1))));//acos(-1));
//   double CEP(acos(-1)-2*acos(-1)*int((energy*origin)/(2*acos(-1))));//acos(-1));

   double pump_probe_delay(8.0/0.02418884);

  double strength_probe(0.000);//0.022);
   double origin_probe(origin+pump_probe_delay);//(4/0.02418884);//3.4/0.02418884);
   double sigma_probe(0.677/0.02418884);//0.887/0.02418884);//0.677/0.02418884);
   double energy_probe(0.05695);//0.08563);//0.05695);
   double CEP_probe(0-2*acos(-1)*int((energy_probe*origin_probe)/(2*acos(-1))));//acos(-1));

   double time(this->m_h*time_index);
   double amplitude(strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*sin(energy*(time-origin)+CEP)
                  +strength_probe*exp(-pow((time-origin_probe),2)/(2*sigma_probe*sigma_probe))*sin(energy_probe*(time-origin_probe)+CEP_probe)
         );
   vector[0]=0;
   vector[1]=0;
   vector[2]=amplitude;
   return 0;

}
