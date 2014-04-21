//@Aimin Yan @5:47 10/14/2005
//After add new point, generate a new pdb file for display

#include "cal_tot_asa_of_res.h"

void regenerate_pdb_file(vector<aa> w)
{    
        int i,j,k;
        k=1;
        ofstream new_pdb("1a04.pdb",ios::out|ios::app);

	for(i=0;i<w.size();i++)
	{
          
     		//cout<<setiosflags(ios::left)<<setw(8)<<w[i].res_num
        	//<<setiosflags(ios::left)<<setw(8)<<w[i].res_id
        	//<<setiosflags(ios::left)<<setw(8)<<w[i].chain_type<<endl;
                    
                for(j=0;j<w[i].atom_str.size();j++)
                {   
     	     new_pdb<<"ATOM  "
                    <<setiosflags(ios::right)<<setw(5)<<k<<"  ";
                
             new_pdb<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(4)<<w[i].atom_str[j].atom_id;

             new_pdb<<setw(4)<<w[i].res_id
                    <<setw(1)<<w[i].chain_type;


             new_pdb<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<w[i].res_num
                    <<setiosflags(ios::fixed)
        	    <<setprecision(3)<<setw(12)<<w[i].atom_str[j].x
        	    <<setprecision(3)<<setw(8)<<w[i].atom_str[j].y
                    <<setprecision(3)<<setw(8)<<w[i].atom_str[j].z
                    <<setprecision(3)<<setw(8)<<w[i].atom_str[j].asa
                    <<setprecision(2)<<setw(6)<<w[i].atom_str[j].b
                    <<endl;
                k=k+1;	
                }

        }//end of for(i=0

       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<10<<setw(5)<<11<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<112<<setw(5)<<113<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<132<<setw(5)<<133<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<185<<setw(5)<<186<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<195<<setw(5)<<196<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<209<<setw(5)<<210<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<248<<setw(5)<<249<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<304<<setw(5)<<305<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<355<<setw(5)<<356<<endl;
       new_pdb<<"CONECT"<<setiosflags(ios::right)<<setw(5)<<444<<setw(5)<<445<<endl;



}

