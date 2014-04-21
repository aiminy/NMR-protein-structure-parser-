//@Aimin Yan
//filename:sum_and_output5.cpp
//output different residue omega value classified by roughness of residue   

#include "cal_tot_asa_of_res.h"

void sum_and_output5(vector<protein> whole_protein)
{
            int i,j,k; 

            char  *temp;

            string residue1[19];
            string temp1;
  
           
            //double angle_total_value,angle_average_value;
            double  angle_total_num,angle_total_num1,angle_total_num2;

            //angle_total_value=0;
            angle_total_num=0;
            angle_total_num1=0;
            angle_total_num2=0;
   
            residue1[0]="ALA";
            residue1[1]="VAL";
            residue1[2]="PHE";
            residue1[3]="PRO";
            residue1[4]="MET";
            residue1[5]="ILE";
            residue1[6]="LEU";
            residue1[7]="ASP";
            residue1[8]="GLU";
            residue1[9]="LYS";
            residue1[10]="ARG";
            residue1[11]="SER";
            residue1[12]="THR";
            residue1[13]="TYR";
            residue1[14]="HIS";
            residue1[15]="CYS";
            residue1[16]="ASN";
            residue1[17]="GLN";
            residue1[18]="TRP";

            temp=new char[20];
                
            double  num_of_aa[19], num_of_aa1[19],num_of_aa2[19];
            double angle_total_value[19],angle_total_value1[19],angle_total_value2[19];
            double angle_average_value[19],angle_average_value1[19],angle_average_value2[19];
            double pro_aa[19],pro_aa1[19],pro_aa2[19];
           
            double t;
            t=0;

            angle_total_num=0;
            angle_total_num1=0;
            angle_total_num2=0;
 
	 for(i=0;i<whole_protein.size();i++)
	 {
           for(j=0;j<whole_protein[i].residue_str.size();j++)
           {
                  //count residue number of 0_1_patch globally 
                if(whole_protein[i].residue_str[j].roughness>0&&whole_protein[i].residue_str[j].roughness<=1)
                  {
                   angle_total_num=angle_total_num+1;
                  }//end if
                  //count residue number of 1_10_patch globally
                if(whole_protein[i].residue_str[j].roughness>1&&whole_protein[i].residue_str[j].roughness<=10)
                  {
                   angle_total_num1=angle_total_num1+1;
                   }//end if

                  //count residue number of  >_10_patch globally  
                if(whole_protein[i].residue_str[j].roughness>10)
                  {
                   angle_total_num2=angle_total_num2+1;
                   }//end if

              }//end for(j=0..
           }//end for(i=0....

        cout<<"global total"<<" "<<angle_total_num<<" "<<angle_total_num1<<" "<<angle_total_num2<<endl;  

       
        for(k=0;k<19;k++)
       {   
          angle_total_value[k]=0; 
          num_of_aa[k]=0;
          angle_average_value[k]=0;

          angle_total_value1[k]=0; 
          num_of_aa1[k]=0;
          angle_average_value1[k]=0;

          angle_total_value2[k]=0; 
          num_of_aa2[k]=0;
          angle_average_value2[k]=0;

          pro_aa[k]=0;
          pro_aa1[k]=0;
          pro_aa2[k]=0;
          

          temp1=residue1[k]+".global_patch_0_1";
          strcpy(temp,temp1.c_str()); 
          //cout<<temp1<<endl;
	  //ofstream file(temp,ios::out|ios::app);
	  //file<<setiosflags(ios::left)<<setw(10)<<"angle"<<"\t"<<"type"<<endl;

	 for(i=0;i<whole_protein.size();i++)
	 {
           for(j=0;j<whole_protein[i].residue_str.size();j++)
           {
                //cout<<whole_protein[i].residue_str[j].res_id<<endl;
	    if(strncmp(whole_protein[i].residue_str[j].res_id,temp,3)==0)
            {
               // cout<<whole_protein[i].residue_str[j].asa<<endl;
               
              if(whole_protein[i].residue_str[j].asa>5)
              {
                //cout<<whole_protein[i].residue_str[j].roughness<<endl;
                   //global center and 0_1_patch
                  
                if(whole_protein[i].residue_str[j].roughness>0&&whole_protein[i].residue_str[j].roughness<=1)
                  {
                   //cout<<"test"<<endl; 
                   angle_total_value[k]=angle_total_value[k]+whole_protein[i].residue_str[j].angle;
                   //angle_total_num=angle_total_num+1;
                   num_of_aa[k]=num_of_aa[k]+1;
                   //file<<setiosflags(ios::left)<<setw(10)<<proteins[i].angle<<"\t"<<proteins[i].res_id<<endl;  
                  }//end if

                  //global center and 1_10_patch
                if(whole_protein[i].residue_str[j].roughness>1&&whole_protein[i].residue_str[j].roughness<=10)
                  {
                   angle_total_value1[k]=angle_total_value1[k]+whole_protein[i].residue_str[j].angle;
                   //angle_total_num1=angle_total_num1+1;
                   num_of_aa1[k]=num_of_aa1[k]+1;
                   //file<<setiosflags(ios::left)<<setw(10)<<proteins[i].angle<<"\t"<<proteins[i].res_id<<endl;  
                   }//end if

                  //global center and >_10_patch  
                if(whole_protein[i].residue_str[j].roughness>10)
                  {
                   angle_total_value2[k]=angle_total_value2[k]+whole_protein[i].residue_str[j].angle;
                   //angle_total_num2=angle_total_num2+1;
                   num_of_aa2[k]=num_of_aa2[k]+1;
                   //file<<setiosflags(ios::left)<<setw(10)<<proteins[i].angle<<"\t"<<proteins[i].res_id<<endl;  
                   }//end if
              }//end if 
  
             }//end if

           }//end for
         }//end for

           angle_average_value[k]=angle_total_value[k]/num_of_aa[k];
           angle_average_value1[k]=angle_total_value1[k]/num_of_aa1[k];
           angle_average_value2[k]=angle_total_value2[k]/num_of_aa2[k];


           //cout<<residue1[k]<<" "<<angle_average_value[k]<<" "<<angle_average_value1[k]<<" "<<angle_average_value2[k]<<endl              
           pro_aa[k]=num_of_aa[k]/angle_total_num*100;
           pro_aa1[k]=num_of_aa1[k]/angle_total_num1*100;
           pro_aa2[k]=num_of_aa2[k]/angle_total_num2*100;
          
           //cout<<residue1[k]<<" "<<num_of_aa[k]<<" "<<num_of_aa1[k]<<" "<<num_of_aa2[k]<<endl;
           //cout<<"patch number"<<angle_total_num<<" "<<angle_total_num1<<" "<<angle_total_num2<<endl;
             
           //cout<<residue1[k]<<" "<<pro_aa[k]<<" "<<pro_aa1[k]<<" "<<pro_aa2[k]<<endl; 
       
      } // end of for k=0
                          
          //ofstream file1_g("aa_g.txt",ios::out);
         // file1_g<<"0_1_patch "<<"\t"<<angle_total_num<<"\t"
           //    <<"1_10_patch"<<"\t"<<angle_total_num1<<"\t"
             //  <<"10_patch  "<<"\t"<<angle_total_num2<<endl;
          
          
         for(i=0;i<19;i++)
          { 
	  cout<<setiosflags(ios::left)<<setw(10)<<residue1[i]<<setiosflags(ios::fixed)
               <<setiosflags(ios::left)<<setprecision(2)<<setw(5)<<pro_aa[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(5)<<pro_aa1[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(5)<<pro_aa2[i]
               <<endl;
          }

          //ofstream file2;
          //file2.open("omega_g.txt",ios::out);
          //file2<<"0_1_patch "<<"\t"<<angle_total_num<<"\t"
           //    <<"1_10_patch"<<"\t"<<angle_total_num1<<"\t"
             //  <<"10_patch  "<<"\t"<<angle_total_num2<<endl;
          
          
         for(i=0;i<19;i++)
          { 
	  cout<<setiosflags(ios::left)<<setw(10)<<residue1[i]<<setiosflags(ios::fixed)
               <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<angle_average_value[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<angle_average_value1[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<angle_average_value2[i]
               <<endl;
          }

          //file1_g.close(); 
          //file2.close(); 

}// end of function sum_and_output5
