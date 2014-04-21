#include "cal_tot_asa_of_res.h"

void output_gene_info_a_protein(vector<protein> whole_protein)     
{    
        int i,j;
	for(i=0;i<whole_protein.size();i++)
	{
          cout<<whole_protein[i].protein_name<<endl;
           cout<<whole_protein[i].residue_str.size()<<endl;

        	for(j=0;j<whole_protein[i].residue_str.size();j++)
        	{
     		cout<<setiosflags(ios::left)<<setw(8)<<whole_protein[i].residue_str[j].res_num
        	<<setiosflags(ios::left)<<setw(8)<<whole_protein[i].residue_str[j].res_id
        	<<setiosflags(ios::left)<<setw(8)<<whole_protein[i].residue_str[j].chain_type
        	<<setiosflags(ios::fixed)
        	<<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[j].angle
        	<<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[j].asa
                <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[j].distance
                <<endl;	
 	        }//end of for

        }//end of for

return;
}

