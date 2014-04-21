#include "cal_tot_asa_of_res.h"

void res_alpha_atom(vector<naa> &w_c3,vector <aa> &w_c2,vector<residue> &a_w_c,vector<aa> w_c)
//void geo_cen_of_res(vector<naa> &w_c3,vector <aa> &w_c2,vector<residue> &a_w_c,vector<aa> w_c)
{
        aa amino;
	amino.res_id       = new char [5];
	
	atom atom_of_res;
	atom_of_res.atom_id = new char [5];
	
	int i,j,k,m;
       
	residue ave_aa;
	
	for(i=0;i<w_c.size();i++)
	 {
	if(strncmp(w_c[i].res_id,"GLY",3)!=0)
	{
	   amino.chain_type=w_c[i].chain_type;
	   amino.res_id=w_c[i].res_id;
	   amino.res_num=w_c[i].res_num;
	   
         //for calpha atom
	   k=0;
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
          if(strncmp(w_c[i].atom_str[j].atom_id,"CA ",3)==0)
	  k=k+1;
	 }//end of for(j=0..
	
	   ave_aa.ave_x=0;
	   ave_aa.ave_y=0;
	   ave_aa.ave_z=0;
	   //if(!amino.atom_str.empty())
	   //{amino.atom_str.clear();}
	   
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
          if(strncmp(w_c[i].atom_str[j].atom_id,"CA ",3)==0)
	   {
	   ave_aa.ave_x=ave_aa.ave_x+w_c[i].atom_str[j].x;
	   ave_aa.ave_y=ave_aa.ave_y+w_c[i].atom_str[j].y;
	   ave_aa.ave_z=ave_aa.ave_z+w_c[i].atom_str[j].z;
	   //atom_of_res.x=atom_of_res.x+w_c[i].atom_str[j].x;
	   //atom_of_res.y=atom_of_res.y+w_c[i].atom_str[j].y;
	   //atom_of_res.z=atom_of_res.z+w_c[i].atom_str[j].z;
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_str[j].asa;
	   //atom_of_res.b=atom_of_res.b+w_c[i].atom_str[j].b;
	   } // end of if
	 } //end of for(j=0..
	
           ave_aa.ave_x=ave_aa.ave_x/k;
	   ave_aa.ave_y=ave_aa.ave_y/k;
	   ave_aa.ave_z=ave_aa.ave_z/k;
	   //ave_aa.asa=ave_aa.asa/k;
	   //ave_aa.ave_b=ave_aa.ave_b/k;
	   //amino.atom_str.push_back(atom_of_res);
 
           a_w_c.push_back(ave_aa); 	 
	   //atom_of_res.atom_id="CA";
           //atom_of_res.x=atom_of_res.x/k;
	   //atom_of_res.y=atom_of_res.y/k;
	   //atom_of_res.z=atom_of_res.z/k;
	   //atom_of_res.asa=atom_of_res.asa/k;
	   //atom_of_res.b=atom_of_res.b/k;
	   //amino.atom_str.push_back(atom_of_res);

         // for c_sc
           
	   //atom_of_res.x=0;
	   //atom_of_res.y=0;
	   //atom_of_res.z=0;
	   //atom_of_res.asa=0;
	   //atom_of_res.b=0;
	  // m=0;
	   //for(j=0;j<w_c[i].atom_str.size();j++)
	 //{
          //if(strncmp(w_c[i].atom_str[j].atom_id,"CA ",3)!=0&&strncmp(w_c[i].atom_str[j].atom_id,"C  ",3)!=0
	  //&&strncmp(w_c[i].atom_str[j].atom_id,"N  ",3)!=0&&strncmp(w_c[i].atom_str[j].atom_id,"O  ",3)!=0)
	  // {
	  // atom_of_res.x=atom_of_res.x+w_c[i].atom_str[j].x;
	  // atom_of_res.y=atom_of_res.y+w_c[i].atom_str[j].y;
	  // atom_of_res.z=atom_of_res.z+w_c[i].atom_str[j].z;
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_str[j].asa;
	  // atom_of_res.b=atom_of_res.b+w_c[i].atom_str[j].b;
	  // m=m+1;
	   //}// end of if(stncmp..
	 //}// end of for(j=0..
	
/*
           if(m!=0) 
	   {atom_of_res.atom_id="C_SC";
	    atom_of_res.x=atom_of_res.x/m;
	    atom_of_res.y=atom_of_res.y/m;
	    atom_of_res.z=atom_of_res.z/m;
	    atom_of_res.asa=atom_of_res.asa/m;
	    atom_of_res.b=atom_of_res.b/m;
	    amino.atom_str.push_back(atom_of_res);
           }// end of if
*/
 
           //if(k!=0)
           //{w_c2.push_back(amino);}
          //end of if(k!=0&&m!=0) 	 
	 } //end of if("!=GLY) 
	 }// end of for(i<w_c.size()..
/*
            
           double c_m_x,c_m_y,c_m_z;
           c_m_x=0;
           c_m_y=0;
           c_m_z=0;

           int toal_aa=0;     

	for(i=0;i<a_w_c.size();i++)
	{
                 c_m_x=c_m_x+a_w_c[i].ave_x;
                 c_m_y=c_m_y+a_w_c[i].ave_y;
                 c_m_z=c_m_z+a_w_c[i].ave_z;
                 total_aa=total+1;
        }//end of for
         
           c_m_x=c_m_x/total_aa;
           c_m_y=c_m_y/total_aa;
           c_m_z=c_m_z/total_aa;
           
           double ave_radius_of_gyration;

         for(i=0;i<a_w_c.size();i++)
         {
          ave_radius_of_gyration=pow((a_w_c[i].ave_x-c_m_x),2)+pow((a_w_c[i].ave_y-c_m_y),2)+pow((a_w_c[i].ave_z),2);
         }//end of for
          
          ave_radius_of_gyration=sqrt(ave_radius_of_gyration/total_aa);
      
          cout<<ave_radius_of_gyration<<endl;
  */
             
}// end of function res_alpha_atom_

