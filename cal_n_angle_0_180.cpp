#include "cal_tot_asa_of_res.h"

void cal_n_angle_0_180(char p_n[],vector<naa>w_c3,char *flag)
{
        
        int    i,n_tot,n_10,n_20,n_30,n_40,n_50,n_60,n_70,n_80,n_90,n_100,s;
        double a_t,a_10,a_20,a_30,a_40,a_50,a_60,a_70,a_80,a_90,a_100;
        double a_av_t,a_av_10,a_av_20,a_av_30,a_av_40,a_av_50,a_av_60,a_av_70,a_av_80,a_av_90,a_av_100;
	      
        vector<std_aa>st_aa_c;	
        char *aa;
           
        produce_std_aa(st_aa_c);

        aa=new char[4];
       for(s=0;s<st_aa_c.size();s++)
	{
         strcpy(aa,st_aa_c[s].res_id);
        n_tot=0;
        n_10=0;
        n_20=0;
        n_30=0;
        n_40=0;
        n_50=0;
        n_60=0;
        n_70=0;
        n_80=0;
        n_90=0;
        n_100=0;
       
        a_t=0;
        a_10=0;
        a_20=0;
        a_30=0;
        a_40=0;
        a_50=0;
        a_60=0;
        a_70=0;
        a_80=0;
        a_90=0;
        a_100=0;
        
        a_av_t=0;
        a_av_10=0;
        a_av_20=0;
        a_av_30=0;
        a_av_40=0;
        a_av_50=0;
        a_av_60=0;
        a_av_70=0;
        a_av_80=0;
        a_av_90=0;
        a_av_100=0;

         ofstream file7(aa,ios::out|ios::app);  
 
        for(i=0;i<w_c3.size();i++)
        {
         if(strncmp(w_c3[i].res_id,aa,3)==0)
	{file7<<setiosflags(ios::left)<<setw(12)<<p_n<<setiosflags(ios::left)<<setw(12)<<
                 w_c3[i].angle<<setiosflags(ios::left)<<setw(12)<<w_c3[i].asa<<endl;
          n_tot=n_tot+1;a_t=a_t+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa<=10.00)
         {n_10=n_10+1;a_10=a_10+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>10.00&&w_c3[i].asa<=20.00)
         {n_20=n_20+1;a_20=a_20+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>20.00&&w_c3[i].asa<=30.00)
         {n_30=n_30+1;a_30=a_30+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>30.00&&w_c3[i].asa<=40.00)
         {n_40=n_40+1;a_40=a_40+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>40.00&&w_c3[i].asa<=50.00)
         {n_50=n_50+1;a_50=a_50+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>50.00&&w_c3[i].asa<=60.00)
         {n_60=n_60+1;a_60=a_60+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>60.00&&w_c3[i].asa<=70.00)
         {n_70=n_70+1;a_70=a_70+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>70.00&&w_c3[i].asa<=80.00)
         {n_80=n_80+1;a_80=a_80+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>80.00&&w_c3[i].asa<=90.00)
         {n_90=n_90+1;a_90=a_90+w_c3[i].angle;}
         
         if(strncmp(w_c3[i].res_id,aa,3)==0&&w_c3[i].asa>90.00)
         {n_100=n_100+1;a_100=a_100+w_c3[i].angle;}
        } // end of for
         
         if(n_tot!=0)
         {a_av_t=a_t/n_tot;}
         
         if(n_10!=0)
         {a_av_10=a_10/n_10;}
         
         if(n_20!=0)
         {a_av_20=a_20/n_20;}
         
         if(n_30!=0)
         {a_av_30=a_30/n_30;}
         
         if(n_40!=0)
         {a_av_40=a_40/n_40;}
         
         if(n_50!=0)
         {a_av_50=a_50/n_50;}
         
         if(n_60!=0)
         {a_av_60=a_60/n_60;}
         
         if(n_70!=0)
         {a_av_70=a_70/n_70;}
         
         if(n_80!=0)
         {a_av_80=a_80/n_80;}
         
         if(n_90!=0)
         {a_av_90=a_90/n_90;}
         
         if(n_100!=0)
         {a_av_100=a_100/n_100;}
        
	ofstream file1,file2;
        char *file1_name;
        char *file2_name;
        file1_name=new char[10];
        file2_name=new char[10];
        
        strcpy(file1_name,strcat(aa,".num"));
        strcpy(file2_name,strtok(aa,"."));
        strcat(file2_name,".angle"); 
            	
	file1.open(file1_name,ios::out|ios::app);
        file2.open(file2_name,ios::out|ios::app);
	
	if(flag)
        {print_head(file1,file2);}
	
	file1<<setiosflags(ios::left)<<setw(12)<<p_n<<setw(8)<<strtok(aa,".")
            <<setiosflags(ios::fixed)
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_tot
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_10
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_20
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_30
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_40
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_50
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_60
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_70
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_80
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_90
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<n_100<<endl;

        file2<<setiosflags(ios::left)<<setw(12)<<p_n<<setw(8)<<strtok(aa,".")
            <<setiosflags(ios::fixed)
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_t
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_10
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_20
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_30
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_40
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_50
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_60
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_70
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_80
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_90
            <<setiosflags(ios::left)<<setprecision(1)<<setw(6)<<a_av_100<<endl;

        file1.close();
        file2.close();   
        file7.close(); 
        aa=new char[4];
       }// end of for(s=0
}
