//	This is linked_list implementation for read file
//	I made four class: Protein,Chain,Residue,Atom
//      for one protein,I made two linked list:
//      atom list:    atom->atom->atom.....atom
//      residue list: residue->residue->residue
//      note in linked_list.h, use linked_list.template, not linked_list.cpp
//      this program read one pdb file, and calculate omega angle, then output *_chain.pdb
//      usage: ./main_read_pdb *.pdb
   
#include "protein.h"

int main(int argc,char *argv[])
{

 if(argc<2)
 {cout<<"usage:"<<"./main_read_pdb asa_file"<<endl;
  exit(0);
 }
 //char* temp;
 //temp=new char;
 //strcpy(temp,argv[1]);
  
 Protein *protein_ptr;
 protein_ptr=new Protein();

 //protein_ptr=new Protein(argv[1],argv[2]);
 protein_ptr=new Protein(argv[1]);
// protein_ptr->SetDifferentAminoAcidTypeInfoForThisProtein();

 //protein_ptr->OutputPrBackBoneAtom();
   //protein_ptr->GetNumberofChainForPr();
  //protein_ptr->SetPrSs(argv[1]);
  //protein_ptr->OutputPrSeSsPo(argv[2]);

     //protein_ptr->OutputPrAaOmegaBasedOnEachChainCetDifferentAminoAcidTypeInfoForThisModel()
     //:wcout<<"test"<<endl;
     //protein_ptr->GetModelPrname();
 //protein_ptr->SetPrBindingSiteAaList();

   /*
   ofstream out_put_file0("small_monomer_p.txt",ios::out|ios::app);
   ofstream out_put_file1("large_monomer_p.txt",ios::out|ios::app);
 
   int num_of_chain,num_of_aa,num_of_atom;
   num_of_chain=protein_ptr->GetNumberofChainForPr();
   num_of_aa=protein_ptr->GetNumberofAaForPr();
   num_of_atom=protein_ptr->GetNumberofAtomForPr();
   
    if(num_of_aa<=100)
    {
     out_put_file0<<temp<<"\t"<<num_of_chain<<"\t"<<num_of_aa<<"\t"<<num_of_atom<<endl;
    }
    else if(num_of_aa>=400)
     {
     out_put_file1<<temp<<"\t"<<num_of_chain<<"\t"<<num_of_aa<<"\t"<<num_of_atom<<endl;
     }
    */

 
 //  protein_ptr->WriteThisProteinSummaryToFile();
  // protein_ptr->CalculateDistanceBetweenCaOfAaForPr();
 //protein_ptr->CalculateDistanceBetweenAtomForPr();

 //cout<<protein_ptr->GetNumberofChainForPr()<<endl;

 /*
 Residue temp_aa;
 Atom temp_atom;
 Protein temp_p;

 temp_aa=protein_ptr->GetPrHeadAa()->data; 

 temp_p=*protein_ptr;


 temp_atom=protein_ptr->GetPrHeadAa()->data.GetAaHeadAtom()->data;

 cout<<"temp_aa"<<temp_aa<<endl; 

 cout<<"temp_atom"<<temp_atom<<endl;

 cout<<"temp_p"<<temp_p<<endl;
*/
 
 // ofstream out_put_file(strcat(strtok(argv[1],"."),".omega_curvature"),ios::out);

/*
  out_put_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<"pr_t\t"
  <<setw(7)<<"pc_x\t"<<setw(7)<<"pc_y\t"<<setw(7)<<"pc_z\t"<<setw(5)<<"aa_t\t"<<setw(5)<<"aa_p\t"
  <<setw(5)<<"aa_c\t"<<setw(7)<<"aa_as\t"<<setw(7)<<"aa_ms\t"<<setw(7)<<"aa_cu\t"<<setw(7)<<"ca_x\t"
  <<setw(7)<<"ca_y"<<setw(7)<<"ca_z"<<setw(7)<<"sc_x"<<setw(7)<<"sc_y"<<setw(7)<<"sc_z\t"
  <<setw(9)<<"aa_omega"<<endl;
 */

  
 //output overall residue 
   //out_put_file<<*protein_ptr;

 //output surface residue
 //protein_ptr->OutputPrSurfaceAalist();

 //output surface convex residue
 //protein_ptr->OutputPrSurfaceConvexAaList();

 //output surface concave residue
 //protein_ptr->OutputPrSurfaceConcaveAaList();

 //output buried residue
 //protein_ptr->OutputPrBurialAaList();

// cout<<"Done "<<protein_ptr->GetPrName()<<endl;
   delete protein_ptr;

}
