//calculate the geometrical center of residue
//then calculate distance between residues

#include "cal_tot_asa_of_res.h"

void geo_cen_of_res(vector<naa> &w_c3,vector <aa> &w_c2,vector<residue> &a_w_c,vector<aa> w_c)
{
        aa amino;
	amino.res_id       = new char [5];
	
	atom atom_of_res;
	atom_of_res.atom_id = new char [5];

        residue ave_aa;
	
	int i,j,k,m;
	
	for(i=0;i<w_c.size();i++)
	 {
	if(strncmp(w_c[i].res_id,"GLY",3)!=0)
	{
	   ave_aa.chain_type=w_c[i].chain_type;
	   ave_aa.res_id=w_c[i].res_id;
	   ave_aa.res_num=w_c[i].res_num;
	   
	   k=0;
           //cout<<w_c[i].atom_str.size()<<endl;
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
	  k=k+1;
	 }
	
	   ave_aa.ave_x=0;
	   ave_aa.ave_y=0;
	   ave_aa.ave_z=0;
	   //ave_aa.ave_b=0;
	   //ave_aa.asa  =0;
	  // if(!amino.atom_str.empty())
	   //{amino.atom_str.clear();}
	   
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
	   ave_aa.ave_x=ave_aa.ave_x+w_c[i].atom_str[j].x;
	   ave_aa.ave_y=ave_aa.ave_y+w_c[i].atom_str[j].y;
	   ave_aa.ave_z=ave_aa.ave_z+w_c[i].atom_str[j].z;
	   //ave_aa.asa=ave_aa.asa+w_c[i].atom_str[j].asa;
           
           //cout<<w_c[i].atom_str[j].b<<endl;
	   //ave_aa.ave_b=ave_aa.ave_b+w_c[i].atom_str[j].b;
	 } //end of for

           //cout<<k<<endl;
           //cout<<ave_aa.ave_x<<endl;
           //cout<<ave_aa.ave_y<<endl;
           //cout<<ave_aa.ave_z<<endl;
           //cout<<ave_aa.asa<<endl;
           //cout<<ave_aa.ave_b<<endl;
            
	   
           ave_aa.ave_x=ave_aa.ave_x/k;
	   ave_aa.ave_y=ave_aa.ave_y/k;
	   ave_aa.ave_z=ave_aa.ave_z/k;
	   //ave_aa.asa=ave_aa.asa/k;
	   //ave_aa.ave_b=ave_aa.ave_b/k;
	   //amino.atom_str.push_back(atom_of_res);
 
           a_w_c.push_back(ave_aa); 	 
	 } //end of if 
	 }// end of for

	for(i=0;i<a_w_c.size();i++)
	{
     		//cout<<setiosflags(ios::left)<<setw(8)<<a_w_c[i].res_num
        	//<<setiosflags(ios::left)<<setw(8)<<a_w_c[i].res_id
        	//<<setiosflags(ios::left)<<setw(8)<<a_w_c[i].chain_type
        	cout<<setiosflags(ios::fixed)
        	<<setiosflags(ios::left)<<setprecision(2)<<setw(7)<<a_w_c[i].ave_x
        	<<setiosflags(ios::left)<<setprecision(2)<<setw(7)<<a_w_c[i].ave_y
        	<<setiosflags(ios::left)<<setprecision(2)<<setw(7)<<a_w_c[i].ave_z
        	//<<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<a_w_c[i].asa
        	//<<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<a_w_c[i].ave_b
                <<endl;	
        }//end of for
 
}
