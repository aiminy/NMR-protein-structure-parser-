#include "cal_tot_asa_of_res.h"

void sum_and_output3(vector<protein> whole_protein)
{
            int i,j,k; 

            //vector< vector<double> > whole_chain;
            //vector<double> exposed;
            //vector<double> interface;
            //vector<double> burial;
            char **residue;
            char *temp;
  
            residue=new char*[20];
           

            for(i=0;i<19;i++)
            {
             residue[i]=new char[4];
            }

            residue[0]="ALA";
            residue[1]="VAL";
            residue[2]="PHE";
            residue[3]="PRO";
            residue[4]="MET";
            residue[5]="ILE";
            residue[6]="LEU";
            residue[7]="ASP";
            residue[8]="GLU";
            residue[9]="LYS";
            residue[10]="ARG";
            residue[11]="SER";
            residue[12]="THR";
            residue[13]="TYR";
            residue[14]="HIS";
            residue[15]="CYS";
            residue[16]="ASN";
            residue[17]="GLN";
            residue[18]="TRP";

            temp=new char[10];
     

       for(k=0;k<19;k++)
   
       {   
          temp=residue[k]; 

	 ofstream file(temp,ios::out|ios::app);
	 file<<setiosflags(ios::left)<<setw(10)<<"angle"<<"\t"<<"type"<<endl;

	 for(i=0;i<whole_protein.size();i++)
	 {
           for(j=0;j<whole_protein[i].residue_str.size();j++)
           {
	    if(strncmp(whole_protein[i].residue_str[j].res_id,temp,3)==0)
            {
              if(whole_protein[i].residue_str[j].asa>5)
              {
               file<<whole_protein[i].residue_str[j].angle<<"\t"<<"E"<<endl; 
              }

             // if(whole_protein[i].residue_str[j].asa<=5&&whole_protein[i].residue_str[j].delta_asa>=1)
             // {
             //  file<<whole_protein[i].residue_str[j].angle<<"\t"<<"I"<<endl; 
             // }
              
	      if(whole_protein[i].residue_str[j].asa<=5)
              {
               file<<whole_protein[i].residue_str[j].angle<<"\t"<<"B"<<endl; 
              }

            }// end of if
           } // end of for
         } // end of for

        file.close();

      } // end of for k=0
                          
}// end of function sum_and_output3 
