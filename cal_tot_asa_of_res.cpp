#include "cal_tot_asa_of_res.h"

void cal_tot_asa_of_res(vector<naa2> &w_c4,vector <aa> w_c)
{
	
	int i,j,k;
        vector <aa_asa> std_asa;
        read_standard_data(std_asa);

        naa2 aa; 
        aa.res_id=new char[5]; 
        double tot_atom_asa;
        
        for(i=0;i<w_c.size();i++)
	{
         aa.res_num=w_c[i].res_num;
         aa.chain_type=w_c[i].chain_type;
         aa.res_id=w_c[i].res_id;

	tot_atom_asa=0;
        for(j=0;j<w_c[i].atom_str.size();j++)
	{
	tot_atom_asa=tot_atom_asa+w_c[i].atom_str[j].asa;
	}// end of for
        
        for(k=0;k<std_asa.size();k++)
        {
         if(strncmp(aa.res_id,std_asa[k].res_id,3)==0)
         {aa.tot_atom_asa=(tot_atom_asa/std_asa[k].asa_1)*100;}
         }//end of if
         w_c4.push_back(aa);
        } //end of for

}
