
#include "cal_tot_asa_of_res.h"

void ter_atom_of_side_chain(vector<naa> &w_c3,vector <aa> &w_c2,vector<aa> w_c,char asa[],char *flag)
{
        aa amino;
	amino.res_id       = new char [5];
	
	atom atom_of_res;
	atom_of_res.atom_id = new char [5];


        vector<std_aa>st_aa_c;	
        char *st_aa,*r_atom;
           
        produce_std_aa(st_aa_c); 
	int i,j,k,m,s;
	
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
	 }
	
	   atom_of_res.x=0;
	   atom_of_res.y=0;
	   atom_of_res.z=0;
	   atom_of_res.asa=0;
	   atom_of_res.b=0;
	   if(!amino.atom_str.empty())
	   {amino.atom_str.clear();}
	   
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
          if(strncmp(w_c[i].atom_str[j].atom_id,"CA ",3)==0)
	   {
	   atom_of_res.x=atom_of_res.x+w_c[i].atom_str[j].x;
	   atom_of_res.y=atom_of_res.y+w_c[i].atom_str[j].y;
	   atom_of_res.z=atom_of_res.z+w_c[i].atom_str[j].z;
	   atom_of_res.asa=atom_of_res.asa+w_c[i].atom_str[j].asa;
	   atom_of_res.b=atom_of_res.b+w_c[i].atom_str[j].b;
	   } // end of if
	 } //end of for
	
	   atom_of_res.atom_id="CA ";
           atom_of_res.x=atom_of_res.x/k;
	   atom_of_res.y=atom_of_res.y/k;
	   atom_of_res.z=atom_of_res.z/k;
	   atom_of_res.asa=atom_of_res.asa/k;
	   atom_of_res.b=atom_of_res.b/k;
	   amino.atom_str.push_back(atom_of_res);

         // for termial atom of side-chain:
         // AA:terminal atom
         // ARG:CZ,LYS:NZ,SER:OG,ILE:CD1,THR:CB,ASN:CG
         // TYR:OH,PRO:CG,LEU:CG,GLU:CD,ASP:CG,ALA:CB
         // CYS:SG,VAL:CB,MET:CE,GLN:CD,HIS:CG,PHE:CZ
         // TRP:CH2
	   
	   atom_of_res.x=0;
	   atom_of_res.y=0;
	   atom_of_res.z=0;
	   atom_of_res.asa=0;
	   atom_of_res.b=0;
	   m=0;
           
           char *pch;
           pch=new char[4];

           st_aa=new char[4];
           r_atom=new char[4];
 
          for(s=0;s<st_aa_c.size();s++)
           {strcpy(st_aa,st_aa_c[s].res_id);
            strcpy(r_atom,st_aa_c[s].r_atom); 
           
	    for(j=0;j<w_c[i].atom_str.size();j++)
	   {
            if(strncmp(w_c[i].res_id,st_aa,3)==0&&strncmp(w_c[i].atom_str[j].atom_id,r_atom,3)==0)
	      {
               atom_of_res.x=atom_of_res.x+w_c[i].atom_str[j].x;
	       atom_of_res.y=atom_of_res.y+w_c[i].atom_str[j].y;
	       atom_of_res.z=atom_of_res.z+w_c[i].atom_str[j].z;
	       atom_of_res.asa=atom_of_res.asa+w_c[i].atom_str[j].asa;
	       atom_of_res.b=atom_of_res.b+w_c[i].atom_str[j].b;
	       m=m+1;
               pch=w_c[i].atom_str[j].atom_id;
	      }// end of if
           }//end of for(j=0..

            st_aa=new char[4];
            r_atom=new char[4];
	
          }//end of for(s=0..
          
          if(m!=0) 
	   {atom_of_res.atom_id=pch;
	    atom_of_res.x=atom_of_res.x/m;
	    atom_of_res.y=atom_of_res.y/m;
	    atom_of_res.z=atom_of_res.z/m;
	    atom_of_res.asa=atom_of_res.asa/m;
	    atom_of_res.b=atom_of_res.b/m;
	    amino.atom_str.push_back(atom_of_res);
           }// end of if
          
         if(k!=0&&m!=0)
	 {w_c2.push_back(amino);}
         	 
	 } //end of if 
	 }// end of for

         double c_x,c_y,c_z;
	 c_x=0;
	 c_y=0;
	 c_z=0;
	 int n;
	 n=0;
	for(i=0;i<w_c2.size();i++)
	{
          for(j=0;j<w_c2[i].atom_str.size();j++)
	 
	{
	c_x=c_x+w_c2[i].atom_str[j].x;
	c_y=c_y+w_c2[i].atom_str[j].y;
	c_z=c_z+w_c2[i].atom_str[j].z;
	n=n+1;
	}// end of for
	
	} //end of for
	c_x=c_x/n;
	c_y=c_y/n;
	c_z=c_z/n;
 
        
        naa amino_n;
        amino_n.res_id=new char[5];
	
        double a_c_x,a_c_y,a_c_z,t1,t2,t3,ta_x,ta_y,ta_z;
	double sc_a_x,sc_a_y,sc_a_z,angle,ac;
        
	a_c_x=0;a_c_y=0;a_c_z=0;t1=0;t2=0;t3=0;
	sc_a_x=0;sc_a_y=0;sc_a_z=0;
        ta_x=0;ta_y=0;ta_z=0;
	
	for(i=0;i<w_c2.size();i++)
	{
	for(j=0;j<w_c2[i].atom_str.size();j++)
	{
 	 if(strncmp(w_c2[i].atom_str[j].atom_id,"CA ",3)==0)
	 {
          t1=sqrt(squ_of_num(w_c2[i].atom_str[j].x-c_x)+squ_of_num(w_c2[i].atom_str[j].y-c_y)+
	          squ_of_num(w_c2[i].atom_str[j].z-c_z));
	  a_c_x=(w_c2[i].atom_str[j].x-c_x)/t1;
	  a_c_y=(w_c2[i].atom_str[j].y-c_y)/t1;
	  a_c_z=(w_c2[i].atom_str[j].z-c_z)/t1;
	  ta_x=w_c2[i].atom_str[j].x;
	  ta_y=w_c2[i].atom_str[j].y;
	  ta_z=w_c2[i].atom_str[j].z;
	  }//end of if  center of structure -> calpha
 	 
	 else{
          t2=sqrt(squ_of_num(w_c2[i].atom_str[j].x-ta_x)+squ_of_num(w_c2[i].atom_str[j].y-ta_y)+
            squ_of_num(w_c2[i].atom_str[j].z-ta_z));
	  sc_a_x=(ta_x-w_c2[i].atom_str[j].x)/t2;
	  sc_a_y=(ta_y-w_c2[i].atom_str[j].y)/t2;
	  sc_a_z=(ta_z-w_c2[i].atom_str[j].z)/t2;
	 }//end of else terminal atom of side chain ->calpha
	}// end of for
          
        t3=sqrt(squ_of_num(a_c_x)+squ_of_num(a_c_y)+squ_of_num(a_c_z))*sqrt(squ_of_num(sc_a_x)+squ_of_num(sc_a_y)+
	        squ_of_num(sc_a_z));
	ac=acos((a_c_x*sc_a_x+a_c_y*sc_a_y+a_c_z*sc_a_z)/t3);
	angle=180*ac/3.14;

	amino_n.res_num=w_c2[i].res_num;
	amino_n.res_id=w_c2[i].res_id;
	amino_n.chain_type=w_c2[i].chain_type;
	amino_n.angle=angle;

	w_c3.push_back(amino_n);
	}//end of for
         
}// end of function geo_of_
