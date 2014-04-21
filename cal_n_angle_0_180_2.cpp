
#include "cal_tot_asa_of_res.h"

void cal_n_angle_0_180_2(char p_n[],vector<naa>w_c3,double a,double b, ofstream &file1, ofstream &file2)
{
         int i;
          
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
	 if(w_c3[i].asa>a&&w_c3[i].asa<=b)
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
         }//end of if
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
	  char buffer[3];
          gcvt(b,3,buffer);
          strcat(buffer,".txt"); 
	  ofstream file(buffer,ios::out|ios::app);
           file<<setiosflags(ios::left)<<setw(11)
	       <<"pn/rn"
               <<setiosflags(ios::left)<<setw(6)
	       <<"ALA"
               <<setiosflags(ios::left)<<setw(6)
	       <<"VAL"
               <<setiosflags(ios::left)<<setw(6)
	       <<"PHE"
               <<setiosflags(ios::left)<<setw(6)
	       <<"PRO"
               <<setiosflags(ios::left)<<setw(6)
	       <<"MET"
               <<setiosflags(ios::left)<<setw(6)
	       <<"ILE"
               <<setiosflags(ios::left)<<setw(6)
	       <<"LEU"
               <<setiosflags(ios::left)<<setw(6)
	       <<"ASP"
               <<setiosflags(ios::left)<<setw(6)
	       <<"GLU"
               <<setiosflags(ios::left)<<setw(6)
	       <<"LYS"
               <<setiosflags(ios::left)<<setw(6)
	       <<"ARG"
               <<setiosflags(ios::left)<<setw(6)
	       <<"SER"
               <<setiosflags(ios::left)<<setw(6)
	       <<"THR"
               <<setiosflags(ios::left)<<setw(6)
	       <<"TYR"
               <<setiosflags(ios::left)<<setw(6)
	       <<"HIS"
               <<setiosflags(ios::left)<<setw(6)
	       <<"CYS"
               <<setiosflags(ios::left)<<setw(6)
	       <<"ASN"
               <<setiosflags(ios::left)<<setw(6)
	       <<"GLN"
               <<setiosflags(ios::left)<<setw(6)
	       <<"TRP"<<endl;
           
	file<<setiosflags(ios::left)<<setw(11)<<p_n
	       <<setiosflags(ios::fixed)
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_ala
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_val
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_phe
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_pro
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_met
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_ile
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_leu
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_asp
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_glu
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_lys
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_arg
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_ser
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_thr
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_tyr
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_his
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_cys
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_asn
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_gln
               <<setiosflags(ios::left)<<setprecision(1)<<setw(6)
	       <<ave_ang_trp<<endl;
             file.close();
}
