#include <vector>


double calculate_mean_of_double_array(vector<double> w);
double calculate_mean_of_double_array(vector<double> w)
{
   double total,average;
   int i,num;
   total=0;average=0;num=0;
   for(i=0;i<w.size();i++)
   {
     total=total+w[i];
     num=num+1;
   }
   average=total/num;
   return average;
}

//a function to calculate standard deviation of a double array 
double calculate_sd_of_double_array(vector<double> w);
double calculate_sd_of_double_array(vector<double> w)
{
 double total,sd,average;
 int i,num;
 total=0;sd=0;average=0;num=0;
 average=calculate_mean_of_double_array(w);
 for(i=0;i<w.size();i++)
 {
  total=total+(w[i]-average)*(w[i]-average);
  num=num+1;
  }
 sd=sqrt(total/(num-1));
 return sd;
}

