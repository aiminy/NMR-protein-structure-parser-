#include "cal_tot_asa_of_res.h"

void read_standard_data(vector<aa_asa> &std_asa)
{
           ifstream infile; 
           infile.open("standard.data");
           aa_asa aa_asa_c; 
           char temp[50]; 
           aa_asa_c.res_id=new char[5];
           double asa1,asa2,asa3,asa4,asa5,asa6,asa7,asa8,asa9,asa10,asa11,asa12,asa13,asa14; 
          
	   int i=0; 
           while(infile.peek()!=-1)
          {
           infile>>temp;
           if(strncmp(temp,"ATOM",4)==0)
           {
           infile.ignore(8);
           
           infile.get(aa_asa_c.res_id,4);
           infile.ignore(2);

           infile.get(temp,7);
           aa_asa_c.asa_1=atof(temp);
           infile.ignore(3);

           infile.get(temp,4);
           aa_asa_c.asa_2=atof(temp);
           infile.ignore(1);
           
           infile.get(temp,7);
           aa_asa_c.asa_3=atof(temp);
           infile.ignore(3);
           
           infile.get(temp,4);
           aa_asa_c.asa_4=atof(temp);
           infile.ignore(1);
           
           infile.get(temp,7);
           aa_asa_c.asa_5=atof(temp);
           infile.ignore(3);
           
           infile.get(temp,4);
           aa_asa_c.asa_6=atof(temp);
           infile.ignore(1);
           
           infile.get(temp,7);
           aa_asa_c.asa_7=atof(temp);
           infile.ignore(3);
           
           infile.get(temp,4);
           aa_asa_c.asa_8=atof(temp);
           infile.ignore(1);
           
           infile.get(temp,7);
           aa_asa_c.asa_9=atof(temp);
           infile.ignore(3);
           
           infile.get(temp,4);
           aa_asa_c.asa_10=atof(temp);
           infile.ignore(1);
           
           infile.get(temp,7);
           aa_asa_c.asa_11=atof(temp);
           infile.ignore(3);
           
           infile.get(temp,4);
           aa_asa_c.asa_12=atof(temp);
           infile.ignore(1);
           
           infile.get(temp,7);
           aa_asa_c.asa_13=atof(temp);
           infile.ignore(3);
           
           infile.get(temp,4);
           aa_asa_c.asa_14=atof(temp);
         
           std_asa.push_back(aa_asa_c);
           aa_asa_c.res_id=new char[5]; 
           }//end of if 
          else{infile.getline(temp,256);}  
          }
          infile.close();

}

