#include "cal_tot_asa_of_res.h"
enum AA {ALA,VAL,PHE,PRO,MET,ILE,LEU,ASP,GLU,LYS,ARG,SER,THR,TYR,HIS,CYS,ASN,GLN,TRP};
int  Radius[]={20,40,60,80,100,120,140,160,180};

void sum_based_on_radius(vector<protein> whole_protein)
{
         int i,j;  
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
         
            
         double angle_total_value=0;
         int angle_total_num=0;

	 for(i=0;i<whole_protein.size();i++)
	 {
           for(j=0;j<whole_protein[i].residue_str.size();j++)
           {  
             // angle_total_value=angle_total_value+whole_protein[i].residue_str[j].angle;
             // angle_total_num=angle_total_num+1;
             // }
             /// ///}
         if(whole_protein[i].residue_str[j].distance>160 
            &&whole_protein[i].residue_str[j].distance<=180) 
         {  
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"ALA",3)==0)
	   {
	   ang_ala=ang_ala+whole_protein[i].residue_str[j].angle;
	   k_ala=k_ala+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"VAL",3)==0)
	   {
	   ang_val=ang_val+whole_protein[i].residue_str[j].angle;
	   k_val=k_val+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"PHE",3)==0)
	   {
	   ang_phe=ang_phe+whole_protein[i].residue_str[j].angle;
	   k_phe=k_phe+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"PRO",3)==0)
	   {
	   ang_pro=ang_pro+whole_protein[i].residue_str[j].angle;
	   k_pro=k_pro+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"MET",3)==0)
	   {
	   ang_met=ang_met+whole_protein[i].residue_str[j].angle;
	   k_met=k_met+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"ILE",3)==0)
	   {
	   ang_ile=ang_ile+whole_protein[i].residue_str[j].angle;
	   k_ile=k_ile+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"LEU",3)==0)
	   {
	   ang_leu=ang_leu+whole_protein[i].residue_str[j].angle;
	   k_leu=k_leu+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"ASP",3)==0)
	   {
	   ang_asp=ang_asp+whole_protein[i].residue_str[j].angle;
	   k_asp=k_asp+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"GLU",3)==0)
	   {
	   ang_glu=ang_glu+whole_protein[i].residue_str[j].angle;
	   k_glu=k_glu+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"LYS",3)==0)
	   {
	   ang_lys=ang_lys+whole_protein[i].residue_str[j].angle;
	   k_lys=k_lys+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"ARG",3)==0)
	   {
	   ang_arg=ang_arg+whole_protein[i].residue_str[j].angle;
	   k_arg=k_arg+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"SER",3)==0)
	   {
	   ang_ser=ang_ser+whole_protein[i].residue_str[j].angle;
	   k_ser=k_ser+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"THR",3)==0)
	   {
	   ang_thr=ang_thr+whole_protein[i].residue_str[j].angle;
	   k_thr=k_thr+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"TYR",3)==0)
	   {
	   ang_tyr=ang_tyr+whole_protein[i].residue_str[j].angle;
	   k_tyr=k_tyr+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"HIS",3)==0)
	   {
	   ang_his=ang_his+whole_protein[i].residue_str[j].angle;
	   k_his=k_his+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"CYS",3)==0)
	   {
	   ang_cys=ang_cys+whole_protein[i].residue_str[j].angle;
	   k_cys=k_cys+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"ASN",3)==0)
	   {
	   ang_asn=ang_asn+whole_protein[i].residue_str[j].angle;
	   k_asn=k_asn+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"GLN",3)==0)
	   {
	   ang_gln=ang_gln+whole_protein[i].residue_str[j].angle;
	   k_gln=k_gln+1;
	   }
	  if(strncmp(whole_protein[i].residue_str[j].res_id,"TRP",3)==0)
	   {
	   ang_trp=ang_trp+whole_protein[i].residue_str[j].angle;
	   k_trp=k_trp+1;
	   }
          }// end of if 
         } //end of for(j=0
         }//end of for(i=0	
          
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
	   
	 ofstream file("radius_160_180.txt",ios::out|ios::app);

           file<<setiosflags(ios::left)<<setw(5)
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
    
           file<<setiosflags(ios::fixed)
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
