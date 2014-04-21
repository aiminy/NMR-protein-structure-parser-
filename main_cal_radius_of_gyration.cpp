//@Aimin Yan
//filename:read_asa.cpp
//usage: ./read_asa *.asa
//read monomer protein file *.asa in this directory calculated from naccess
//calculate omega angle and indentify SASA of residue and atom
//output to file angle_asa.txt
//angle_asa.txt format:
//res_num res_name chain_id res_omega res_rel_asa
//for example:
//  1       GLU       A       67.00     34.87 
 
#include "cal_tot_asa_of_res.h"

int main (int argc, char *argv[]) {
	
	coord_struct  asa;
	int    i,j;
	char p_n[20];
	vector<char> chain;
        vector<aa> w_c, w_c2;	
        vector<aa_asa> std_asa; 
        vector<naa2>w_c4;
        vector<naa>w_c3;

        vector<residue>a_w_c;
	
	//vector<complex_aa> total_aa;   
	
	char *c_a;
	c_a=new char[1];
	

         vector<protein> whole_protein;
         vector<aa> wc_c_aa2,wc_t_aa2;
         vector<naa> wc_c_aa3,wc_t_aa3;
	
        
	 ifstream protein_file("monomer.txt",ios::in);

         char *file_name;
         file_name=new char[10];

         while(protein_file>>file_name)
        {
         //cout<<setiosflags(ios::left)<<setw(10)<<file_name
           //  <<setiosflags(ios::left)<<setw(10)<<"Rg is"; 

	//read *.asa file and produce struct asa,vector w_c 
	read_asa_file(w_c,asa,file_name);
        cal_tot_asa_of_res(w_c4,w_c);

        // cout<<setiosflags(ios::left)<<setw(15)<<file_name<<setiosflags(ios::left)<<setw(10)<<w_c.size()<<endl;
        // w_c.clear();

        //geo_cen_of_res(wc_c_aa3,wc_c_aa2,a_w_c,w_c);
          //res_alpha_atom(wc_c_aa3,wc_c_aa2,a_w_c,w_c);
        
        //cout<<whole_protein.size()<<endl; 
        cal_radius_of_gyration(w_c3,w_c2,a_w_c,w_c);


            
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
         for(j=0;j<a_w_c.size();j++)
         {
         if(wc_c_aa3[i].res_num==a_w_c[j].res_num&&wc_c_aa3[i].chain_type==a_w_c[j].chain_type&&
            (strncmp(wc_c_aa3[i].res_id,a_w_c[j].res_id,3)==0))
          {
           wc_c_aa3[i].distance=a_w_c[j].distance;
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
          one_residue.distance=wc_c_aa3[i].distance;
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
        a_w_c.clear();
        }//end of while 
         //end of read different protein 
	

        //cout<<whole_protein.size()<<endl; 

          
         output_gene_info_a_protein(whole_protein); 
         sum_based_on_radius(whole_protein);
         sum_and_output(whole_protein);
         sum_and_output3(whole_protein); 
        // sum_count_output(whole_protein); 
       
         protein_file.close(); 
      
     
        return(0);
     }//end of main
