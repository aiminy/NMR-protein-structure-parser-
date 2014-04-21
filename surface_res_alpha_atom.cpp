#include "cal_tot_asa_of_res.h"

void surface_res_alpha_atom(vector<naa> &w_c3,vector <aa> &w_c2,vector<aa> w_c,char asa[],char *flag)
{
        aa amino;
	amino.res_id       = new char [5];
	
	atom atom_of_res;
	atom_of_res.atom_id = new char [5];
	
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
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
          if(strncmp(w_c[i].atom_str[j].atom_id,"CA ",3)==0)
	  k=k+1;
	 }//end of for(j=0..
	
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
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_str[j].asa;
	   //atom_of_res.b=atom_of_res.b+w_c[i].atom_str[j].b;
	   } // end of if
	 } //end of for(j=0..
	
	   atom_of_res.atom_id="CA";
           atom_of_res.x=atom_of_res.x/k;
	   atom_of_res.y=atom_of_res.y/k;
	   atom_of_res.z=atom_of_res.z/k;
	   //atom_of_res.asa=atom_of_res.asa/k;
	   //atom_of_res.b=atom_of_res.b/k;
	   amino.atom_str.push_back(atom_of_res);

         // for c_sc
           
	   //atom_of_res.x=0;
	   //atom_of_res.y=0;
	   //atom_of_res.z=0;
	   atom_of_res.asa=0;
	   //atom_of_res.b=0;
	   m=0;
	   for(j=0;j<w_c[i].atom_str.size();j++)
	 {
          //if(strncmp(w_c[i].atom_str[j].atom_id,"CA ",3)!=0&&strncmp(w_c[i].atom_str[j].atom_id,"C  ",3)!=0
	  //&&strncmp(w_c[i].atom_str[j].atom_id,"N  ",3)!=0&&strncmp(w_c[i].atom_str[j].atom_id,"O  ",3)!=0)
	  // {
	  // atom_of_res.x=atom_of_res.x+w_c[i].atom_str[j].x;
	  // atom_of_res.y=atom_of_res.y+w_c[i].atom_str[j].y;
	  // atom_of_res.z=atom_of_res.z+w_c[i].atom_str[j].z;
	   atom_of_res.asa=atom_of_res.asa+w_c[i].atom_str[j].asa;
	  // atom_of_res.b=atom_of_res.b+w_c[i].atom_str[j].b;
	  // m=m+1;
	   }// end of if(stncmp..
	 }// end of for(j=0..
	
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
 
           if(k!=0&&m!=0)
           {w_c2.push_back(amino);}
          //end of if(k!=0&&m!=0) 	 
	 } //end of if("!=GLY) 
	 }// end of for(i<w_c.size()..

         
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
 	 if(strncmp(w_c2[i].atom_str[j].atom_id,"CA",2)==0)
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
 	 
	if(strncmp(w_c2[i].atom_str[j].atom_id,"C_SC",4)==0)
	 {
          t2=sqrt(squ_of_num(w_c2[i].atom_str[j].x-ta_x)+squ_of_num(w_c2[i].atom_str[j].y-ta_y)+
            squ_of_num(w_c2[i].atom_str[j].z-ta_z));
	  sc_a_x=(ta_x-w_c2[i].atom_str[j].x)/t2;
	  sc_a_y=(ta_y-w_c2[i].atom_str[j].y)/t2;
	  sc_a_z=(ta_z-w_c2[i].atom_str[j].z)/t2;
	 }//end of if    center of side chain ->calpha
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

         double ang_ala,ave_ang_ala,k_ala;
         double ang_val,ave_ang_val,k_val;
         double ang_phe,ave_ang_phe,k_phe;
         double ang_pro,ave_ang_pro,k_pro;
         double ang_met,ave_ang_met,k_met;
         double ang_ile,ave_ang_ile,k_ile;
         double ang_leu,ave_ang_leu,k_leu;
         double ang_asp,ave_ang_asp,k_asp;
         double ang_glu,ave_ang_glu,k_glu;
         double ang_lys,ave_ang_lys,k_lys;
         double ang_arg,ave_ang_arg,k_arg;
         double ang_ser,ave_ang_ser,k_ser;
         double ang_thr,ave_ang_thr,k_thr;
         double ang_tyr,ave_ang_tyr,k_tyr;
         double ang_his,ave_ang_his,k_his;
         double ang_cys,ave_ang_cys,k_cys;
         double ang_asn,ave_ang_asn,k_asn;
         double ang_gln,ave_ang_gln,k_gln;
         double ang_trp,ave_ang_trp,k_trp;

	 ang_ala=0;k_ala=0;ave_ang_ala=0;
	 ang_val=0;k_val=0;ave_ang_val=0;
	 ang_phe=0;k_phe=0;ave_ang_phe=0;
	 ang_pro=0;k_pro=0;ave_ang_pro=0;
	 ang_met=0;k_met=0;ave_ang_met=0;
	 ang_ile=0;k_ile=0;ave_ang_ile=0;
	 ang_leu=0;k_leu=0;ave_ang_leu=0;
	 ang_asp=0;k_asp=0;ave_ang_asp=0;
	 ang_glu=0;k_glu=0;ave_ang_glu=0;
	 ang_lys=0;k_lys=0;ave_ang_lys=0;
	 ang_arg=0;k_arg=0;ave_ang_arg=0;
	 ang_ser=0;k_ser=0;ave_ang_ser=0;
	 ang_thr=0;k_thr=0;ave_ang_thr=0;
	 ang_tyr=0;k_tyr=0;ave_ang_tyr=0;
	 ang_his=0;k_his=0;ave_ang_his=0;
	 ang_cys=0;k_cys=0;ave_ang_cys=0;
	 ang_asn=0;k_asn=0;ave_ang_asn=0;
	 ang_gln=0;k_gln=0;ave_ang_gln=0;
	 ang_trp=0;k_trp=0;ave_ang_trp=0;
         
	 for(i=0;i<w_c3.size();i++)
	 {
	  if(strncmp(w_c3[i].res_id,"ALA",3)==0)
	   {
	   ang_ala=ang_ala+w_c3[i].angle;
	   k_ala=k_ala+1;
	   }
	  if(strncmp(w_c3[i].res_id,"VAL",3)==0)
	   {
	   ang_val=ang_val+w_c3[i].angle;
	   k_val=k_val+1;
	   }
	  if(strncmp(w_c3[i].res_id,"PHE",3)==0)
	   {
	   ang_phe=ang_phe+w_c3[i].angle;
	   k_phe=k_phe+1;
	   }
	  if(strncmp(w_c3[i].res_id,"PRO",3)==0)
	   {
	   ang_pro=ang_pro+w_c3[i].angle;
	   k_pro=k_pro+1;
	   }
	  if(strncmp(w_c3[i].res_id,"MET",3)==0)
	   {
	   ang_met=ang_met+w_c3[i].angle;
	   k_met=k_met+1;
	   }
	  if(strncmp(w_c3[i].res_id,"ILE",3)==0)
	   {
	   ang_ile=ang_ile+w_c3[i].angle;
	   k_ile=k_ile+1;
	   }
	  if(strncmp(w_c3[i].res_id,"LEU",3)==0)
	   {
	   ang_leu=ang_leu+w_c3[i].angle;
	   k_leu=k_leu+1;
	   }
	  if(strncmp(w_c3[i].res_id,"ASP",3)==0)
	   {
	   ang_asp=ang_asp+w_c3[i].angle;
	   k_asp=k_asp+1;
	   }
	  if(strncmp(w_c3[i].res_id,"GLU",3)==0)
	   {
	   ang_glu=ang_glu+w_c3[i].angle;
	   k_glu=k_glu+1;
	   }
	  if(strncmp(w_c3[i].res_id,"LYS",3)==0)
	   {
	   ang_lys=ang_lys+w_c3[i].angle;
	   k_lys=k_lys+1;
	   }
	  if(strncmp(w_c3[i].res_id,"ARG",3)==0)
	   {
	   ang_arg=ang_arg+w_c3[i].angle;
	   k_arg=k_arg+1;
	   }
	  if(strncmp(w_c3[i].res_id,"SER",3)==0)
	   {
	   ang_ser=ang_ser+w_c3[i].angle;
	   k_ser=k_ser+1;
	   }
	  if(strncmp(w_c3[i].res_id,"THR",3)==0)
	   {
	   ang_thr=ang_thr+w_c3[i].angle;
	   k_thr=k_thr+1;
	   }
	  if(strncmp(w_c3[i].res_id,"TYR",3)==0)
	   {
	   ang_tyr=ang_tyr+w_c3[i].angle;
	   k_tyr=k_tyr+1;
	   }
	  if(strncmp(w_c3[i].res_id,"HIS",3)==0)
	   {
	   ang_his=ang_his+w_c3[i].angle;
	   k_his=k_his+1;
	   }
	  if(strncmp(w_c3[i].res_id,"CYS",3)==0)
	   {
	   ang_cys=ang_cys+w_c3[i].angle;
	   k_cys=k_cys+1;
	   }
	  if(strncmp(w_c3[i].res_id,"ASN",3)==0)
	   {
	   ang_asn=ang_asn+w_c3[i].angle;
	   k_asn=k_asn+1;
	   }
	  if(strncmp(w_c3[i].res_id,"GLN",3)==0)
	   {
	   ang_gln=ang_gln+w_c3[i].angle;
	   k_gln=k_gln+1;
	   }
	  if(strncmp(w_c3[i].res_id,"TRP",3)==0)
	   {
	   ang_trp=ang_trp+w_c3[i].angle;
	   k_trp=k_trp+1;
	   }
         }//end of for	
          
	   if(k_ala!=0)
	   {ave_ang_ala=ang_ala/k_ala;}
	   if(k_val!=0)
	   {ave_ang_val=ang_val/k_val;}
	   if(k_phe!=0)
	   {ave_ang_phe=ang_phe/k_phe;}
	   if(k_pro!=0)
	   {ave_ang_pro=ang_pro/k_pro;}
	   if(k_met!=0)
	   {ave_ang_met=ang_met/k_met;}
	   if(k_ile!=0)
	   {ave_ang_ile=ang_ile/k_ile;}
	   if(k_leu!=0)
	   {ave_ang_leu=ang_leu/k_leu;}
	   if(k_asp!=0)
	   {ave_ang_asp=ang_asp/k_asp;}
	   if(k_glu!=0)
	   {ave_ang_glu=ang_glu/k_glu;}
	   if(k_lys!=0)
	   {ave_ang_lys=ang_lys/k_lys;}
	   if(k_arg!=0)
	   {ave_ang_arg=ang_arg/k_arg;}
	   if(k_ser!=0)
	   {ave_ang_ser=ang_ser/k_ser;}
	   if(k_thr!=0)
	   {ave_ang_thr=ang_thr/k_thr;}
	   if(k_tyr!=0)
	   {ave_ang_tyr=ang_tyr/k_tyr;}
	   if(k_his!=0)
	   {ave_ang_his=ang_his/k_his;}
	   if(k_cys!=0)
	   {ave_ang_cys=ang_cys/k_cys;}
	   if(k_asn!=0)
	   {ave_ang_asn=ang_asn/k_asn;}
	   if(k_gln!=0)
	   {ave_ang_gln=ang_gln/k_gln;}
	   if(k_trp!=0)
	   {ave_ang_trp=ang_trp/k_trp;}
	   
	  ofstream file("angle_based_center.txt",ios::out|ios::app);

	if(flag)
         {
           file<<setiosflags(ios::left)<<setw(13)
	       <<"pn/rn"
               <<setiosflags(ios::left)<<setw(5)
	       <<"ALA"
               <<setiosflags(ios::left)<<setw(5)
	       <<"VAL"
               <<setiosflags(ios::left)<<setw(5)
	       <<"PHE"
               <<setiosflags(ios::left)<<setw(5)
	       <<"PRO"
               <<setiosflags(ios::left)<<setw(5)
	       <<"MET"
               <<setiosflags(ios::left)<<setw(5)
	       <<"ILE"
               <<setiosflags(ios::left)<<setw(5)
	       <<"LEU"
               <<setiosflags(ios::left)<<setw(5)
	       <<"ASP"
               <<setiosflags(ios::left)<<setw(5)
	       <<"GLU"
               <<setiosflags(ios::left)<<setw(5)
	       <<"LYS"
               <<setiosflags(ios::left)<<setw(5)
	       <<"ARG"
               <<setiosflags(ios::left)<<setw(5)
	       <<"SER"
               <<setiosflags(ios::left)<<setw(5)
	       <<"THR"
               <<setiosflags(ios::left)<<setw(5)
	       <<"TYR"
               <<setiosflags(ios::left)<<setw(5)
	       <<"HIS"
               <<setiosflags(ios::left)<<setw(5)
	       <<"CYS"
               <<setiosflags(ios::left)<<setw(5)
	       <<"ASN"
               <<setiosflags(ios::left)<<setw(5)
	       <<"GLN"
               <<setiosflags(ios::left)<<setw(5)
	       <<"TRP"<<endl;
	    }
    
           file<<setiosflags(ios::left)<<setw(13)<<asa
	       <<setiosflags(ios::fixed)
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_ala
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_val
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_phe
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_pro
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_met
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_ile
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_leu
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_asp
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_glu
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_lys
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_arg
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_ser
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_thr
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_tyr
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_his
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_cys
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_asn
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_gln
               <<setiosflags(ios::left)<<setprecision(0)<<setw(5)
	       <<ave_ang_trp<<endl;
             file.close();
}// end of function geo_of_

