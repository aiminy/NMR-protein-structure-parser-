//@Aimin Yan
//8/23/2005
//filename:read_asa_main.cpp
//usage: ./read_asa_main monomer.ruf.pstrm
//output to file results3_asa.txt
 
#include "cal_tot_asa_of_res.h"

int main (int argc, char *argv[]) {

          
	if (argc < PDB_INDEX+1) 
	{
	cout <<"Usage:"<<argv[0]<<" monomer.ruf.pstrm\n";
	exit(1);
	}//end of if
	
	coord_struct  asa;
	int    i,j;
	char p_n[20];
	vector<char> chain;
        vector<aa> w_c, w_c2;	

        vector<aa_asa> std_asa; 
        vector<naa2>w_c4;
        vector<naa>w_c3;

        vector<residue>a_w_c;
	
        vector<protein_rou> proteins_rou;	
	
	char *c_a;
	c_a=new char[1];
	

         vector<protein> whole_protein;
         vector<aa> wc_c_aa2,wc_t_aa2;
         vector<naa> wc_c_aa3,wc_t_aa3;
	
        
	 ifstream protein_file("monomer_u.txt",ios::in);

         char *file_name;
         file_name=new char[10];

         while(protein_file>>file_name)
        {

	read_asa_file(w_c,asa,file_name);
        cal_tot_asa_of_res(w_c4,w_c);

        // cout<<setiosflags(ios::left)<<setw(15)<<file_name<<setiosflags(ios::left)<<setw(10)<<w_c.size()<<endl;
        // w_c.clear();

        //geo_cen_of_res(wc_c_aa3,wc_c_aa2,a_w_c,w_c);
          //res_alpha_atom(wc_c_aa3,wc_c_aa2,a_w_c,w_c);
        
        //cout<<whole_protein.size()<<endl; 
           

            
	geo_cen_of_side_chain2(wc_c_aa3,wc_c_aa2,w_c);
	ter_atom_of_side_chain2(wc_t_aa3,wc_t_aa2,w_c);
 

	for(i=0;i<wc_c_aa3.size();i++)
        {
           wc_c_aa3[i].angle_t=1000;
	}

	
	for(i=0;i<wc_c_aa3.size();i++)
        {
         for(j=0;j<wc_t_aa3.size();j++)
         {
         if(wc_c_aa3[i].res_num==wc_t_aa3[j].res_num&&wc_c_aa3[i].chain_type==wc_t_aa3[j].chain_type&&
            (strncmp(wc_c_aa3[i].res_id,wc_t_aa3[j].res_id,3)==0))
          {
           wc_c_aa3[i].angle_t=wc_t_aa3[j].angle;
          }//end of if
         }//end of for
        }// end of for
        
	for(i=0;i<wc_c_aa3.size();i++)
        {
         for(j=0;j<w_c4.size();j++)
         {
         if(wc_c_aa3[i].res_num==w_c4[j].res_num&&wc_c_aa3[i].chain_type==w_c4[j].chain_type&&
            (strncmp(wc_c_aa3[i].res_id,w_c4[j].res_id,3)==0))
          {
           wc_c_aa3[i].asa=w_c4[j].tot_atom_asa;
          }//end of if
         }//end of for
        }// end of for
/*
	cout<<setiosflags(ios::left)<<setw(8)<<"res_num"
        <<setiosflags(ios::left)<<setw(8)<<"res_id"
        <<setiosflags(ios::left)<<setw(8)<<"ch_type"
        <<setiosflags(ios::left)<<setw(8)<<"an_c_c"
        <<setiosflags(ios::left)<<setw(8)<<"an_c_t"
	<<setiosflags(ios::left)<<setw(8)<<"asa"<<endl;
	
	for(i=0;i<wc_c_aa3.size();i++)
	{
	cout<<setiosflags(ios::left)<<setw(8)<<wc_c_aa3[i].res_num
        <<setiosflags(ios::left)<<setw(8)<<wc_c_aa3[i].res_id
        <<setiosflags(ios::left)<<setw(8)<<wc_c_aa3[i].chain_type
        <<setiosflags(ios::fixed)
        <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<wc_c_aa3[i].angle
        <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<wc_c_aa3[i].angle_t
 	<<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<wc_c_aa3[i].asa<<endl;
	}//end of for
  */
        
        
        residue one_residue;
        one_residue.res_id=new char[5];

        
        protein one_protein;
        one_protein.protein_name=new char[10];

        strcpy(one_protein.protein_name,strtok(file_name,"."));

        for(i=0;i<wc_c_aa3.size();i++)
        {
          one_residue.res_id=wc_c_aa3[i].res_id;
          one_residue.res_num=wc_c_aa3[i].res_num;
          one_residue.chain_type=wc_c_aa3[i].chain_type;
          one_residue.angle=wc_c_aa3[i].angle;
          one_residue.asa=wc_c_aa3[i].asa;
          one_protein.residue_str.push_back(one_residue); 
        }
        whole_protein.push_back(one_protein);

        one_protein.residue_str.clear();
        wc_c_aa2.clear();
        wc_t_aa2.clear();
        wc_c_aa3.clear();
        wc_t_aa3.clear();
        w_c4.clear();
        w_c.clear();  
         
        }//end of while 

        //these two functions are in read_dimer.cpp file     
	read_pdb_file(proteins_rou,argv[PDB_INDEX]);
        calculate_roughness_of_residue(proteins_rou);   
	
	 ofstream output("results3_asa.txt",ios::out);
         int m,n,tt;

	for(i=0;i<whole_protein.size();i++)
         {
          for(m=0;m<whole_protein[i].residue_str.size();m++)
           {whole_protein[i].residue_str[m].roughness=-100;}
         }

          for(j=0;j<proteins_rou.size();j++)
            {
             for(n=0;n<proteins_rou[j].w_c_rou.size();n++)
              {
              proteins_rou[j].w_c_rou[n].asa=0;  
              }
            }               
       

         tt=0;    
	for(i=0;i<whole_protein.size();i++)
	{
          for(j=0;j<proteins_rou.size();j++)
           {
            if(strcmp(whole_protein[i].protein_name,strtok(proteins_rou[j].protein_rou_id,"."))==0) 
           {
        	for(m=0;m<whole_protein[i].residue_str.size();m++)
        	{
                  for(n=0;n<proteins_rou[j].w_c_rou.size();n++)
                   {
                     if(whole_protein[i].residue_str[m].res_num==proteins_rou[j].w_c_rou[n].res_num&&
                       *whole_protein[i].residue_str[m].res_id==*proteins_rou[j].w_c_rou[n].res_id&&
                        whole_protein[i].residue_str[m].chain_type==*proteins_rou[j].w_c_rou[n].chain_type&&
                        whole_protein[i].residue_str[m].asa>5&&
                        proteins_rou[j].w_c_rou[n].roughness>0
                       )
                     {

                     output<<setiosflags(ios::left)<<setw(10)<<whole_protein[i].protein_name
                           <<setiosflags(ios::left)<<setw(6)<<whole_protein[i].residue_str[m].res_num
                           <<setiosflags(ios::left)<<setw(4)<<whole_protein[i].residue_str[m].res_id
                           <<setiosflags(ios::left)<<setw(3)<<whole_protein[i].residue_str[m].chain_type
                           <<setiosflags(ios::left)<<setw(10)<<proteins_rou[j].protein_rou_id
                           <<setiosflags(ios::left)<<setw(6)<<proteins_rou[j].w_c_rou[n].res_num
                           <<setiosflags(ios::left)<<setw(4)<<proteins_rou[j].w_c_rou[n].res_id
                           <<setiosflags(ios::left)<<setw(3)<<proteins_rou[j].w_c_rou[n].chain_type
                           <<setiosflags(ios::fixed)
        	     <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[m].angle
        	     <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[m].asa
                     <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<proteins_rou[j].w_c_rou[n].roughness<<endl;
                     whole_protein[i].residue_str[m].roughness=proteins_rou[j].w_c_rou[n].roughness;
                     proteins_rou[j].w_c_rou[n].asa=whole_protein[i].residue_str[m].asa;
                     tt=tt+1; 
                     }//end if
                   }//end for
 	        }//end of for
            }//end of if
           }//end of for 
        }//end of for
	 
         ofstream output1("results3_asa1.txt",ios::out);

	for(i=0;i<whole_protein.size();i++)
         {
          for(m=0;m<whole_protein[i].residue_str.size();m++)
           {
            if(whole_protein[i].residue_str[m].roughness>0)
            {
                     output1<<setiosflags(ios::left)<<setw(10)<<whole_protein[i].protein_name
                           <<setiosflags(ios::left)<<setw(6)<<whole_protein[i].residue_str[m].res_num
                           <<setiosflags(ios::left)<<setw(4)<<whole_protein[i].residue_str[m].res_id
                           <<setiosflags(ios::left)<<setw(3)<<whole_protein[i].residue_str[m].chain_type
                           <<setiosflags(ios::fixed)
        	     <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[m].angle
        	     <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<whole_protein[i].residue_str[m].asa
                     <<whole_protein[i].residue_str[m].roughness<<endl;
            } 
           }
         }
 
         ofstream output2("results3_asa2.txt",ios::out);

          for(j=0;j<proteins_rou.size();j++)
            {
             for(n=0;n<proteins_rou[j].w_c_rou.size();n++)
              {
               if(proteins_rou[j].w_c_rou[n].asa>5)
                 {                
                     output2<<setiosflags(ios::left)<<setw(10)<<proteins_rou[j].protein_rou_id
                           <<setiosflags(ios::left)<<setw(6)<<proteins_rou[j].w_c_rou[n].res_num
                           <<setiosflags(ios::left)<<setw(4)<<proteins_rou[j].w_c_rou[n].res_id
                           <<setiosflags(ios::left)<<setw(3)<<proteins_rou[j].w_c_rou[n].chain_type
                           <<setiosflags(ios::fixed)
                           <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<proteins_rou[j].w_c_rou[n].roughness
                           <<setiosflags(ios::left)<<setprecision(2)<<setw(8)<<proteins_rou[j].w_c_rou[n].asa<<endl;
                 }         
              }
            }               

         sum_and_output(whole_protein);
         sum_and_output3(whole_protein); 
         sum_and_output5(whole_protein); 
      
         cout<<"total surface residue(asa>5 and roughness>0) is:"<<tt<<endl;
         protein_file.close(); 
      
     
        return(0);
     }//end of main
