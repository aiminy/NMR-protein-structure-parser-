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
#include "linked_list.h"

int main (int argc, char *argv[]) {

 
        if(argc<2)
        {cout<<"usage:main_calculate_distance file"<<endl;
         exit(1);
        }
	

        Protein pro;
         
        node<Protein> *head_protein,*previous_protein;
 
        node<Atom> *head_atom;  
	
        node<Residue> *head_aa;  
        
	 ifstream protein_file(argv[1],ios::in);

         char *file_name;
         file_name=new char[10];

         head_aa=NULL;
         head_protein=NULL;
         previous_protein=NULL;
 
         while(protein_file>>file_name)
        {
        Protein *pro;
        pro=new Protein;
        read_based_atom(pro,file_name);  
        insert_an_node(head_protein,previous_protein,pro);
        

        }//end of while 

        
         node<Protein>* cursor_protein;
         node<Residue>* cursor_aa;
         node<atom>*    cursor_atom;
 
         cursor_protein=head_protein;

         while(cursor_protein!=NULL)
         {
           cout<<cursor_protein->data.GetPrName()<<endl;             
           cursor_aa=cursor_protein->data.GetPrHeadAa();

           inter_atom_distance_matrix2(cursor_aa);
/*
           while(cursor_aa!=NULL)
           {
             if(cursor_aa->data.GetAaCalphaAtom()->data.GetAtom_asa()>5.0)
             { 
               cout<<cursor_aa->data.GetAaCalphaAtom()->data<<endl;
             } //end if
           cursor_aa=cursor_aa->next; 
           } //end while(cursor_aa
*/

           cursor_protein=cursor_protein->next;  
          }//end while(cursor_protein
        


         protein_file.close(); 
                

        return(0);
     }//end of main
