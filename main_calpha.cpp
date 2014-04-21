//@Aimin Yan
//filename:main_calpha.cpp
//to compile:g++ -o main_calpha main_calpha.cpp
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

	if (argc < PDB_INDEX+1) 
	{
	cout <<"Usage: prompt> "<<argv[0]<<" file.pdb\n";
	exit(1);
	}//end of if
	
        vector<protein_rou> proteins_rou;	
             
        read_pdb_file(proteins_rou,argv[PDB_INDEX]);
        calculate_roughness_of_residue(proteins_rou);
        cout<<"ok"<<endl;
       
     
        return(0);
     }//end of main
