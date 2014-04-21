//@Aimin Yan @5:29 10/15/2005
//calculate whole structure center: ws
//calculate center of triangle made by three closest calpha atom to each surface calpha atom(asa>5%):tc
//calculate normal vector to above triangle:tc_n;
//calculate vecotr w_s_t: ws->tc
//calculate angle between w_s_t and tc_n;
//output this angle

#include "cal_tot_asa_of_res.h"

void output_global_angle(vector<aa> w)
{    
        int i;
        ofstream new_pdb("1a04.global",ios::out|ios::app);

	for(i=0;i<w.size();i++)
	{
          if(w[i].global_angle>0)
           {
     	     new_pdb<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<w[i].res_num
                    <<setw(4)<<w[i].res_id
                    <<setw(1)<<w[i].chain_type;

             new_pdb<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)
                    <<setiosflags(ios::fixed)
        	    <<setprecision(3)<<setw(12)<<w[i].global_angle
                    <<endl;
           }

        }//end of for(i=0

}

