/*
this program calculate the center of each surface patch
then based on this center of the patch to calculate omega
angle 
*/

#include "cal_tot_asa_of_res.h"

void geo_cen_of_side_chain3(vector<aa_rou> &w_c3,vector<aa_rou> &w_c2,vector<aa_rou> w_c);

void geo_cen_of_side_chain3(vector<aa_rou> w_c3,vector <aa_rou> w_c2,vector<aa_rou> w_c)
{
       
        aa_rou amino;
	amino.res_id       = new char [5];
	
	atom_rou atom_of_res;
	atom_of_res.atom_rou_id = new char [5];
	
	int i,j,k,m;
	
	for(i=0;i<w_c.size();i++)
	{
	    if(strncmp(w_c[i].res_id,"GLY",3)!=0)
	   {
	    amino.chain_type=w_c[i].chain_type;
	    amino.res_id=w_c[i].res_id;
	    amino.res_num=w_c[i].res_num;
	   
           //for calpha atom
	   k=0;
	   for(j=0;j<w_c[i].atom_rou_str.size();j++)
	   {
            if(strncmp(w_c[i].atom_rou_str[j].atom_rou_id,"CA ",3)==0)
	    k=k+1;
	   }//end of for(j=0..
	
	   atom_of_res.x=0;
	   atom_of_res.y=0;
	   atom_of_res.z=0;
	   //atom_of_res.asa=0;
	   atom_of_res.b=0;
	   if(!amino.atom_rou_str.empty())
	   {amino.atom_rou_str.clear();}
	   
	   for(j=0;j<w_c[i].atom_rou_str.size();j++)
	 {
          if(strncmp(w_c[i].atom_rou_str[j].atom_rou_id,"CA ",3)==0)
	   {
	   atom_of_res.x=atom_of_res.x+w_c[i].atom_rou_str[j].x;
	   atom_of_res.y=atom_of_res.y+w_c[i].atom_rou_str[j].y;
	   atom_of_res.z=atom_of_res.z+w_c[i].atom_rou_str[j].z;
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_rou_str[j].asa;
	   atom_of_res.b=atom_of_res.b+w_c[i].atom_rou_str[j].b;
	   } // end of if
	 } //end of for(j=0..
	
	   atom_of_res.atom_rou_id="CA";
           atom_of_res.x=atom_of_res.x/k;
	   atom_of_res.y=atom_of_res.y/k;
	   atom_of_res.z=atom_of_res.z/k;
	   //atom_of_res.asa=atom_of_res.asa/k;
	   atom_of_res.b=atom_of_res.b/k;
	   amino.atom_rou_str.push_back(atom_of_res);

         // for c_sc
           
	   atom_of_res.x=0;
	   atom_of_res.y=0;
	   atom_of_res.z=0;
	   //atom_of_res.asa=0;
	   atom_of_res.b=0;
	   m=0;
	   for(j=0;j<w_c[i].atom_rou_str.size();j++)
	  {
          if(strncmp(w_c[i].atom_rou_str[j].atom_rou_id,"CA ",3)!=0&&strncmp(w_c[i].atom_rou_str[j].atom_rou_id,"C  ",3)!=0
	  &&strncmp(w_c[i].atom_rou_str[j].atom_rou_id,"N  ",3)!=0&&strncmp(w_c[i].atom_rou_str[j].atom_rou_id,"O  ",3)!=0)
	   {
	   atom_of_res.x=atom_of_res.x+w_c[i].atom_rou_str[j].x;
	   atom_of_res.y=atom_of_res.y+w_c[i].atom_rou_str[j].y;
	   atom_of_res.z=atom_of_res.z+w_c[i].atom_rou_str[j].z;
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_rou_str[j].asa;
	   atom_of_res.b=atom_of_res.b+w_c[i].atom_rou_str[j].b;
	   m=m+1;
	   }// end of if(stncmp..
	  }// end of for(j=0..
	
           if(m!=0) 
	   {atom_of_res.atom_rou_id="C_SC";
	    atom_of_res.x=atom_of_res.x/m;
	    atom_of_res.y=atom_of_res.y/m;
	    atom_of_res.z=atom_of_res.z/m;
	    //atom_of_res.asa=atom_of_res.asa/m;
	    atom_of_res.b=atom_of_res.b/m;
	    amino.atom_rou_str.push_back(atom_of_res);
           }// end of if
 
           if(k!=0&&m!=0)
           {w_c2.push_back(amino);}//end of if(k!=0&&m!=0) 	 
	 }//end of if("!=GLY)
	
        }// end of for(i<w_c.size()..

         
         double c_x,c_y,c_z;
	 c_x=0;
	 c_y=0;
	 c_z=0;
	 int n;
	 n=0;
	for(i=0;i<w_c2.size();i++)
	{
          for(j=0;j<w_c2[i].atom_rou_str.size();j++)
	 
	{
	c_x=c_x+w_c2[i].atom_rou_str[j].x;
	c_y=c_y+w_c2[i].atom_rou_str[j].y;
	c_z=c_z+w_c2[i].atom_rou_str[j].z;
	n=n+1;
	}// end of for(j=0;
	
	} //end of for(i=0;
	c_x=c_x/n;
	c_y=c_y/n;
	c_z=c_z/n;
 
        
        aa_rou amino_n;
        amino_n.res_id=new char[5];
	
        double a_c_x,a_c_y,a_c_z,t1,t2,t3,ta_x,ta_y,ta_z;
	double sc_a_x,sc_a_y,sc_a_z,angle,ac;
        
	a_c_x=0;a_c_y=0;a_c_z=0;t1=0;t2=0;t3=0;
	sc_a_x=0;sc_a_y=0;sc_a_z=0;
        ta_x=0;ta_y=0;ta_z=0;
	
	for(i=0;i<w_c2.size();i++)
	{
	for(j=0;j<w_c2[i].atom_rou_str.size();j++)
	{
 	 if(strncmp(w_c2[i].atom_rou_str[j].atom_rou_id,"CA",2)==0)
	 {
          t1=sqrt(squ_of_num(w_c2[i].atom_rou_str[j].x-c_x)+squ_of_num(w_c2[i].atom_rou_str[j].y-c_y)+
	          squ_of_num(w_c2[i].atom_rou_str[j].z-c_z));
	  a_c_x=(w_c2[i].atom_rou_str[j].x-c_x)/t1;
	  a_c_y=(w_c2[i].atom_rou_str[j].y-c_y)/t1;
	  a_c_z=(w_c2[i].atom_rou_str[j].z-c_z)/t1;
	  ta_x=w_c2[i].atom_rou_str[j].x;
	  ta_y=w_c2[i].atom_rou_str[j].y;
	  ta_z=w_c2[i].atom_rou_str[j].z;
	  }//end of if  center of structure -> calpha
 	 
	if(strncmp(w_c2[i].atom_rou_str[j].atom_rou_id,"C_SC",4)==0)
	 {
          t2=sqrt(squ_of_num(w_c2[i].atom_rou_str[j].x-ta_x)+squ_of_num(w_c2[i].atom_rou_str[j].y-ta_y)+
            squ_of_num(w_c2[i].atom_rou_str[j].z-ta_z));
	  sc_a_x=(ta_x-w_c2[i].atom_rou_str[j].x)/t2;
	  sc_a_y=(ta_y-w_c2[i].atom_rou_str[j].y)/t2;
	  sc_a_z=(ta_z-w_c2[i].atom_rou_str[j].z)/t2;
	 }//end of if    center of side chain ->calpha
	}// end of for(j=0;j<w_c2[i]
          
        t3=sqrt(squ_of_num(a_c_x)+squ_of_num(a_c_y)+squ_of_num(a_c_z))*sqrt(squ_of_num(sc_a_x)+squ_of_num(sc_a_y)+
	        squ_of_num(sc_a_z));
	ac=acos((a_c_x*sc_a_x+a_c_y*sc_a_y+a_c_z*sc_a_z)/t3);
	angle=180*ac/3.14;

	amino_n.res_num=w_c2[i].res_num;
	amino_n.res_id=w_c2[i].res_id;
	amino_n.chain_type=w_c2[i].chain_type;
	amino_n.angle=angle;

	w_c3.push_back(amino_n);
	}//end of for(i=0;i<w_c2.....)

}// end of function geo_of_
