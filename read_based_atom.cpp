//	This is linked_list implementation for read ASA file
//	I made a Atom class, insert each atom into a linked list
//	example:atom->atom->.....atom
// 	this implementation is better than array or vector implementation
//	11:01PM 10/24/2005 @ Aimin Yan
//      note: do make clean first, then do make, then function change, otherwise keep old function

#include "aa.h"
#include "protein.h"
//#include "linked_list.h"
#include <fstream>

using namespace std;

#define Record_Max 90
//template <class Item>
void read_based_atom(Protein*& protein,char *file)
{

 ifstream input(file,ios::in);

 if (input.bad()) 
 {
  cerr<<"Error: Unable to open "<<file<<endl; 
  exit (2);
 }//end of if

 char *temp;
 temp=new char;
 
 Atom atom;
 Residue *aa;
 aa=new Residue;

 //Protein *pro;
 //pro=new Protein;

  protein->SetPrName(file);

 
 node<Atom> *aa_head_atom,*aa_tail_atom;
 aa_head_atom=NULL;aa_tail_atom=NULL;

 node<Residue> *head_aa,*previous_aa;
 head_aa=NULL;
 previous_aa=NULL;

 //node<Protein> *previous_pro;
 //previous_pro=NULL;

 
 node<Atom> *previous_atom;
 previous_atom=NULL;
 
 while(input>>temp)                  //read into ATOM ,but not white space 
  {
   if(strncmp(temp,"ATOM",4)!=0)
   {input.getline(temp,Record_Max);} //if not start from ATOM, read away whold line, include newline(\n)
   else
   {
    input.ignore(2);                 //throw away two white space ~~
    input.get(temp,6);               //read 5 character, addtion 1 for null character(\0) ~~~~1\0
    atom.SetAtomIndex(atoi(temp));
    input.ignore(2);                 //~~
    input.get(temp,4);               //N~~
    atom.SetAtomName(temp);
    input.ignore(1);                 //~
    input.get(temp,4);               //GLU\0
    atom.SetResidueName(temp);
    //aa.SetAaName(temp);
    input.ignore(1);                 //~
    input.get(temp,2);               //A\0
    atom.SetChainType(temp);
    //aa.SetAaChain(temp);
    input.get(temp,5);               //~~~5\0
    atom.SetResidueNum(atoi(temp));
    //aa.SetAaNum(atoi(temp));
    input.ignore(4);                 //~~~~
    input.get(temp,9);               //~~28.823\0
    atom.SetAtom_x(atof(temp));
    input.get(temp,9);               //~-10.606
    atom.SetAtom_y(atof(temp));
    input.get(temp,9);               //~~39.660\0
    atom.SetAtom_z(atof(temp));
    input.get(temp,9);               //~~44.853\0
    atom.SetAtom_asa(atof(temp));
    input.get(temp,7);               //~~1.65\0
    atom.SetAtom_radius(atof(temp));

    //insert_an_node(head_atom,previous_atom,atom);
    //cout<<"test2"<<atom.GetResidueName()<<endl;   
   
    //cout<<"atom1"<<endl;
    //cout<<atom<<endl;
   
     
    if(aa->GetAaHeadAtom()==NULL||
       strncmp(aa->GetAaTailAtom()->data.GetResidueName(),atom.GetResidueName(),3)==0&&
       strncmp(aa->GetAaTailAtom()->data.GetChainType(),atom.GetChainType(),2)==0    &&
       aa->GetAaTailAtom()->data.GetResidueNum()==atom.GetResidueNum())
    {
     aa->SetAaAtomList(atom);
    }
    else if(head_aa==NULL||previous_aa->data.GetAaNum()!=aa->GetAaNum())
    {

     insert_an_node(head_aa,previous_aa,aa);
     aa=new Residue;
     aa->SetAaAtomList(atom);
     //cout<<*aa;
     
    }
   }//end else 
  }//end while
     insert_an_node(head_aa,previous_aa,aa);
   

   //display_node_list(aa_head_atom);
   //display_node_list(head_aa);

     node<Residue> *cursor;
     cursor=head_aa;
   
     while(cursor!=NULL)
     {
     protein->SetPrAaList(cursor->data);
     cursor=cursor->next;
     }

     //cout<<*protein<<endl;
   
    //insert_an_node(head_protein,previous_pro,pro);
        
    //cout<<*pro<<endl;
/*
      cout<<head_aa->data.GetAaNum()<<endl;
      cout<<previous_aa->data.GetAaNum()<<endl;
      cout<<aa->GetAaNum()<<endl;      
 */
     
    //cout<<*aa;

    //cout<<aa->GetAaHeadAtom()->data.GetResidueNum()<<endl;

}
