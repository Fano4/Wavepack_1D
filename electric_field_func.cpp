#include <cmath>


double hamilton_matrix::electric_field(double time_index,double* vector)
{
   double strength(this->m_pump_strength);
   double origin(this->m_pump_origin);
   double sigma(this->m_pump_sigma);
   double energy(this->m_pump_energy);
   double CEP(this->m_pump_CEP-2*acos(-1)*int((this->m_pump_energy*this->m_pump_origin)/(2*acos(-1))));

   double pump_probe_delay(this->m_pprobe_delay);

   double strength_probe(this->m_probe_strength);
   double origin_probe(this->m_pump_origin+this->m_pprobe_delay);
   double sigma_probe(this->m_probe_sigma);
   double energy_probe(this->m_probe_energy);
   double CEP_probe(this->m_probe_CEP-2*acos(-1)*int((this->m_probe_energy*origin_probe)/(2*acos(-1))));

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
   double strength(this->m_pump_strength);
   double origin(this->m_pump_origin);
   double sigma(this->m_pump_sigma);
   double energy(this->m_pump_energy);
   double CEP(this->m_pump_CEP-2*acos(-1)*int((energy*origin)/(2*acos(-1))));

   double pump_probe_delay(this->m_pprobe_delay);

   double strength_probe(this->m_probe_strength);
   double origin_probe(origin+pump_probe_delay);
   double sigma_probe(this->m_probe_sigma);
   double energy_probe(this->m_probe_energy);
   double CEP_probe(this->m_probe_CEP-2*acos(-1)*int((energy_probe*origin_probe)/(2*acos(-1))));

   double time(this->m_h*time_index);
   double amplitude(-strength*exp(-pow((time-origin),2)/(2*sigma*sigma))*sin(energy*(time-origin)+CEP)
                  -strength_probe*exp(-pow((time-origin_probe),2)/(2*sigma_probe*sigma_probe))*sin(energy_probe*(time-origin_probe)+CEP_probe)
         );
   vector[0]=0;
   vector[1]=0;
   vector[2]=amplitude;
   return 0;

}
