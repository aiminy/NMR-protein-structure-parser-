#include "model.h"
#include "aa_summary.h"

#define Record_Max 80

Model::Model()
{
 p_name=NULL;
 model_name=NULL;
 model_no=-1;
 num_of_aa=0;
 num_of_atom=0;
 head_atom_ptr=NULL;
 head_aa_ptr=NULL;
 tail_aa_ptr=NULL;
 Center_Of_Structure=NULL;
 //head_chain_ptr=NULL;
 //tail_chain_ptr=NULL;

}

Model::Model(char *Input_File)
{
  /*
 ifstream input(Input_File,ios::in);
 p_name=strtok(Input_File,".");

 char *temp;
 temp=new char[81];
 Atom atom;
 node<Atom> *head_atom,*previous_atom;
 previous_atom=NULL;

 
 while(input>>temp)                  //read into ATOM ,but not white space 
  {
   if(strncmp(temp,"ATOM",4)==0)
   {
    atom.SetAtomType(temp);
   // cout<<temp<<endl;
    input.ignore(3);                 //throw away one white space ~~
    
    input.get(temp,5);               //read 4 character, addtion 1 for null character(\0) ~~~1\0
    atom.SetAtomIndex(atoi(temp));
   // cout<<temp<<endl;
    
    input.ignore(2);                 //~~
    input.get(temp,4);               //N~~
    atom.SetAtomName(temp);
   // cout<<temp<<endl;
    
    input.ignore(1);                 //~
    input.get(temp,4);               //GLU\0
    atom.SetResidueName(temp);
   // cout<<temp<<endl;
    
    input.get();                 //~
    input.get(temp,2);                //read chain
    atom.SetChainType(temp);
   // cout<<temp<<endl;
    input.get(); 
 
    input.get(temp,4);               //246\0
    atom.SetResidueNum(atoi(temp));
   // cout<<temp<<endl;
    //input.get(temp,5);               //~~~5\0
    //atom.SetResidueNum(atoi(temp));
    input.ignore(5);                 //~~~~
    input.get(temp,8);               //~~~2.980\0
    atom.SetAtom_x(atof(temp));
   // cout<<temp<<endl;
    
    input.ignore(1);
    input.get(temp,8);               //~14.360\0
    atom.SetAtom_y(atof(temp));
   // cout<<temp<<endl;
    
    input.ignore(1);
    input.get(temp,8);               //~~20.984\0
    atom.SetAtom_z(atof(temp));
   // cout<<temp<<endl;
    //input.get(temp,9);               //~~20.984\0
    //atom.SetAtom_asa(atof(temp));
    input.ignore(2);
    input.get(temp,5);
    atom.SetAtomOcc(atof(temp));
    //cout<<temp<<endl;
    
    //input.ignore(1);
    input.get(temp,7);               //~~~33.53\0
    atom.SetAtomBfactor(atof(temp));
    //cout<<temp<<endl;
    
    input.ignore(6);
    input.get(temp,5);               //~~14.99\0
    atom.SetProteinName(temp);
    //cout<<temp<<endl;
    
    input.get(temp,5);               //~-0.35\0
    atom.SetLineNumber(atoi(temp));
    //cout<<temp<<endl;
    insert_an_node(head_atom,previous_atom,atom);
   }//end else
   else if (strncmp(temp,"HETATM",6)==0)
    { 
    atom.SetAtomType(temp);
    //cout<<temp<<endl;
    input.ignore(1);                 //throw away one white space ~~
    
    input.get(temp,5);               //read 4 character, addtion 1 for null character(\0) ~~~1\0
    atom.SetAtomIndex(atoi(temp));
    //cout<<temp<<endl;
    
    input.ignore(2);                 //~~
    input.get(temp,4);               //N~~
    atom.SetAtomName(temp);
    //cout<<temp<<endl;
    
    input.ignore(1);                 //~
    input.get(temp,4);               //GLU\0
    atom.SetResidueName(temp);
    //cout<<temp<<endl;
    
    //input.ignore(3);                 //~
    input.get();                 //~
    input.get(temp,2);                //read chain
    atom.SetModelType(temp);
    input.get(); 

    input.get(temp,4);               //246\0
    atom.SetResidueNum(atoi(temp));
    //cout<<temp<<endl;
    //input.get(temp,5);               //~~~5\0
    //atom.SetResidueNum(atoi(temp));
    input.ignore(5);                 //~~~~
    input.get(temp,8);               //~~~2.980\0
    atom.SetAtom_x(atof(temp));
   // cout<<temp<<endl;
    
    input.ignore(1);
    input.get(temp,8);               //~14.360\0
    atom.SetAtom_y(atof(temp));
   //cout<<temp<<endl;
    
    input.ignore(1);
    input.get(temp,8);               //~~20.984\0
    atom.SetAtom_z(atof(temp));
    //cout<<temp<<endl;
    //input.get(temp,9);               //~~20.984\0
    //atom.SetAtom_asa(atof(temp));
    input.ignore(2);
    input.get(temp,5);
    atom.SetAtomOcc(atof(temp));
    //cout<<temp<<endl;
    
    //input.ignore(1);
    input.get(temp,7);               //~~~33.53\0
    atom.SetAtomBfactor(atof(temp));
    //cout<<temp<<endl;
    
    input.ignore(6);
    input.get(temp,5);               //~~14.99\0
    atom.SetProteinName(temp);
    //cout<<temp<<endl;
    
    input.get(temp,5);               //~-0.35\0
    atom.SetLineNumber(atoi(temp));
    //cout<<temp<<endl;
    insert_an_node(head_atom,previous_atom,atom);
   }
   else
   {
      input.get(temp,Record_Max);
   }

    //insert_an_node(head_atom,previous_atom,atom);
   

  }//end while
   //display_node_list(head_atom);

  
   head_atom_ptr=head_atom;


   //create residue list of this protein: example: residue_1->residue_2->....
    node<Residue> *cursor_aa_ptr;
    cursor_aa_ptr=head_aa;
    Residue *chain;
    chain=new Model; 
    chain->SetModelAaList(cursor_aa_ptr->data);
    node<Model> *head_chain,*previous_chain;
    previous_chain=NULL;
    cursor_aa_ptr=cursor_aa_ptr->next;

    while(cursor_atom_ptr!=NULL)
   {
     if(aa->GetAaNum()==cursor_atom_ptr->data.GetResidueNum()&&
        strcmp(aa->GetAaName(),cursor_atom_ptr->data.GetResidueName())==0&&
        strcmp(aa->GetAaModel(),cursor_atom_ptr->data.GetModelType())==0)
          {
           aa->SetAaAtomList(cursor_atom_ptr->data);
          }
     else
       {
        insert_an_node(head_chain,previous_chain,*chain);
        chain=new Model;
        chain->SetModelAaList(cursor_aa_ptr->data);
       }
       cursor_aa_ptr=cursor_aa_ptr->next;
    }
        insert_an_node(head_chain,previous_chain,*chain);
        head_chain_ptr=head_chain;
        
       // CalculateOmegaAngleForAaOfThisModel();
*/
}

char* Model::GetModelPrname()
{
return p_name;
}

char* Model::GetModelName()
{
 //return head_aa_ptr->data.GetAaModel();
 return model_name;
}

int Model::GetModelNum()
{
 //return head_aa_ptr->data.GetAaModel();
 return model_no;
}

int Model::GetModelAaNum()
{
 return num_of_aa;
}

int Model::GetModelAtomNum()
{
 return num_of_atom;
}

node<Atom>* Model::GetModelCenter()
{
  Center_Of_Structure=new node<Atom>;
  int AtomId,AaId;
  char *AtomName,*AaName,*Model;

  AtomName=new char;
  AaName=new char;
  Model=new char;
  
  double n,av_x,av_y,av_z,av_occ,av_B,av_asa,av_radius;
  
  n=0;av_x=0;av_y=0;av_z=0;av_asa=0;av_radius=0;
  
  node<Residue> *cursor_aa_ptr;
  cursor_aa_ptr=head_aa_ptr;

  node<Atom> *cursor_atom_ptr;

  while(cursor_aa_ptr!=NULL)
  {
   cursor_atom_ptr=cursor_aa_ptr->data.GetAaHeadAtom();
   while(cursor_atom_ptr!=NULL)
   {
        //only consider ATOM type
     if(strncmp(cursor_atom_ptr->data.GetAtomType(),"ATOM",4)==0)
      { 
       av_x=av_x+cursor_atom_ptr->data.GetAtom_x();
       av_y=av_y+cursor_atom_ptr->data.GetAtom_y();
       av_z=av_z+cursor_atom_ptr->data.GetAtom_z();
       n=n+1;
      }
    cursor_atom_ptr=cursor_atom_ptr->next; 
   }
    cursor_aa_ptr=cursor_aa_ptr->next; 
 }

    av_x=av_x/n;
    av_y=av_y/n;
    av_z=av_z/n;

     
    Center_Of_Structure->data.SetAtomIndex(0);
    Center_Of_Structure->data.SetAtomName("C_C");
    Center_Of_Structure->data.SetResidueName("C_C");
    Center_Of_Structure->data.SetChainType("C_C");
    Center_Of_Structure->data.SetResidueNum(0);    
    Center_Of_Structure->data.SetAtom_x(av_x);
    Center_Of_Structure->data.SetAtom_y(av_y);
    Center_Of_Structure->data.SetAtom_z(av_z);
    Center_Of_Structure->data.SetAtom_asa(0);
    Center_Of_Structure->data.SetAtomRadius(0);

    return Center_Of_Structure;
}

node<Residue>* Model::GetModelHeadAa()
{
 return head_aa_ptr;
}

node<Residue>* Model::GetModelTailAa()
{
 return tail_aa_ptr;
}

node<Atom>* Model::GetModelHeadAtom()
{
return head_aa_ptr->data.GetAaHeadAtom();
}

node<Atom>* Model::GetModelTailAtom()
{
return tail_aa_ptr->data.GetAaTailAtom();
}

node<Atom>* Model::GetModelSurfaceHeadCalphaAtom()
{

}

double Model::GetModelAaOmegaAverage(vector<double> w )
{
   double total,average;
   int i,num;
   
   total=0;average=0;num=0;

   for(i=0;i<w.size();i++)
   {
     total=total+w[i];
     num=num+1;
   }

   if(num!=0)
   {average=total/num;}

   return average;

}

double Model::GetModelAaOmegaSd_Of_Ave(vector<double> w)
{
 double total,sd,average;
 int i,num;
 
 total=0;sd=0;average=0;num=0;

 average=GetModelAaOmegaAverage(w);

 for(i=0;i<w.size();i++)
 {
  total=total+(w[i]-average)*(w[i]-average);
  num=num+1;
  }

 if(num!=0&&num!=1)
 {sd=sqrt(total/(num-1));}

 return sd;

}

int Model::GetNumberofModelForPr()
{

 //int num_of_chain;
 //node<Residue> *cursor_aa_ptr;
 //cursor_aa_ptr=head_aa_ptr;
 //node<char> *cursor_chain_ptr;

 //insert_an_node(head_chain_ptr,tail_chain_ptr,*cursor_aa_ptr->data.GetAaModel());
 
 // while(cursor_aa_ptr!=NULL)
 // {
 // if(*cursor_aa_ptr->data.GetAaModel()!=tail_chain_ptr->data)
 //   {insert_an_node(head_chain_ptr,tail_chain_ptr,*cursor_aa_ptr->data.GetAaModel());}
 //    cursor_aa_ptr=cursor_aa_ptr->next;
  //}

 //num_of_chain=linked_list_length(head_chain_ptr);
// display_node_list(head_chain_ptr);
 

 //cout<<num_of_chain<<",";

 //cursor_chain_ptr=head_chain_ptr;
 //while(cursor_chain_ptr!=NULL)
 //{
 // if(cursor_chain_ptr->next!=NULL)
 // {
  //cout<<cursor_chain_ptr->data<<",";
  //}
  //else
  //{cout<<cursor_chain_ptr->data;}
 //cursor_chain_ptr=cursor_chain_ptr->next;
 //}
 //cout<<endl;
// return num_of_chain;

}


node<Aa_Summary>* Model::GetHeadAa_Summary()
{
return  aa_sum;

}
  

int Model::GetNumberofAaForModel()
{
  int num_of_aa;
  num_of_aa=linked_list_length(head_aa_ptr);
  return num_of_aa;

}
  
int Model::GetNumberofAtomForModel()
{
  int num_of_atom;
  num_of_atom=linked_list_length(head_atom_ptr);
  return num_of_atom;
}

void Model::SetPrName(char *input)
{
  p_name=new char;
  strcpy(p_name,input);
}

void Model::SetModelName(char *input)
{
  model_name=new char;
  strcpy(model_name,input);
}

void Model::SetModelNum(int num)
{
 model_no=num;
}


void Model::SetModelAaNum(int num)
{
  num_of_aa=num;
} 

void Model::SetModelAtomNum(int num)
{
  num_of_atom=num;
}

void Model::SetModelAaList(node<Residue> *residue)
{
 //insert_an_node(head_aa_ptr,tail_aa_ptr,residue);
 head_aa_ptr=residue;
}

void Model::SetModelCenter()
{

}

void Model::SetModelAtomList(node<Atom> *source_atom)
{
 head_atom_ptr=source_atom; 
}

//void Model::SetPrAaList(node<Residue> *source_aa)
//{
// head_aa_ptr=source_aa; 
//}


void Model::SetModelSurfaceAaList()
{
    node<Residue> *head_aa,*previous_aa,*cursor_aa_ptr;
    previous_aa=NULL;
    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(cursor_aa_ptr->data.GetAaAs()>=5)
     {
     insert_an_node(head_aa,previous_aa,cursor_aa_ptr->data);
     }
     cursor_aa_ptr=cursor_aa_ptr->next;
    }
     surface_head_aa_ptr=head_aa;
}

void Model::SetModelSurfaceConvexAaList()
{
    node<Residue> *head_aa,*previous_aa,*cursor_aa_ptr;
    previous_aa=NULL;
    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(cursor_aa_ptr->data.GetAaAs()>=5&&
        cursor_aa_ptr->data.GetAaAverageCurvature()<0)
     {
     insert_an_node(head_aa,previous_aa,cursor_aa_ptr->data);
     }
     cursor_aa_ptr=cursor_aa_ptr->next;
    }
     surface_Convex_head_aa_ptr=head_aa;
}

void Model::SetModelSurfaceConcaveAaList()
{
    node<Residue> *head_aa,*previous_aa,*cursor_aa_ptr;
    previous_aa=NULL;
    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(cursor_aa_ptr->data.GetAaAs()>=5&&
        cursor_aa_ptr->data.GetAaAverageCurvature()>=0)
     {
     insert_an_node(head_aa,previous_aa,cursor_aa_ptr->data);
     }
     cursor_aa_ptr=cursor_aa_ptr->next;
    }
     surface_Concave_head_aa_ptr=head_aa;
}

void Model::SetModelBurialAaList()
{
    node<Residue> *head_aa,*previous_aa,*cursor_aa_ptr;
    previous_aa=NULL;
    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(cursor_aa_ptr->data.GetAaAs()<5)
     {
     insert_an_node(head_aa,previous_aa,cursor_aa_ptr->data);
     }
     cursor_aa_ptr=cursor_aa_ptr->next;
    }
     burial_head_aa_ptr=head_aa;

}

void Model::SetModelConvexAaList()
{
    node<Residue> *head_aa,*previous_aa,*cursor_aa_ptr;
    previous_aa=NULL;
    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(cursor_aa_ptr->data.GetAaAverageCurvature()<0)
     {insert_an_node(head_aa,previous_aa,cursor_aa_ptr->data);}
     cursor_aa_ptr=cursor_aa_ptr->next;
    }
     convex_head_aa_ptr=head_aa;

}

void Model::SetModelConcaveAaList()
{
    node<Residue> *head_aa,*previous_aa,*cursor_aa_ptr;
    previous_aa=NULL;
    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(cursor_aa_ptr->data.GetAaAverageCurvature()>=0)
  {insert_an_node(head_aa,previous_aa,cursor_aa_ptr->data);}
     cursor_aa_ptr=cursor_aa_ptr->next;
    }
     concave_head_aa_ptr=head_aa;

}

void Model::SetModelBindingSiteAaList()
{
  
  node<Atom> *cursor_atom_ptr1;
  node<Atom> *cursor_atom_ptr2;
  node<Atom> *cursor_atom_ptr;
  node<Atom> *head_atom,*previous_atom;
  previous_atom=NULL;

  double distance=0;

  cursor_atom_ptr1=head_atom_ptr;
  while(cursor_atom_ptr1!=NULL)
{    
   if(strncmp(cursor_atom_ptr1->data.GetAtomType(),"ATOM",4)==0)
   {
   cursor_atom_ptr2=head_atom_ptr;
  while(cursor_atom_ptr2!=NULL)
     {
      if(strncmp(cursor_atom_ptr2->data.GetAtomType(),"HETATM",6)==0
                              &&strncmp(cursor_atom_ptr2->data.GetAtomName(),"HOH",3)!=0)
     {
     distance=sqrt(pow((cursor_atom_ptr1->data.GetAtom_x()-
                        cursor_atom_ptr2->data.GetAtom_x()),2)+
                   pow((cursor_atom_ptr1->data.GetAtom_y()-
                        cursor_atom_ptr2->data.GetAtom_y()),2)+
                   pow((cursor_atom_ptr1->data.GetAtom_z()-
                        cursor_atom_ptr2->data.GetAtom_z()),2));
     if(distance<4.0&&previous_atom==NULL)
       {insert_an_node(head_atom,previous_atom,cursor_atom_ptr1->data);}
     else if(distance<4.0&&strcmp(previous_atom->data.GetResidueName(),cursor_atom_ptr1->data.GetResidueName())!=0
                         &&previous_atom->data.GetResidueNum()!=cursor_atom_ptr1->data.GetResidueNum())
       {insert_an_node(head_atom,previous_atom,cursor_atom_ptr1->data);}
     }
     cursor_atom_ptr2=cursor_atom_ptr2->next;
   }// end of 2rd while
    }  
     cursor_atom_ptr1=cursor_atom_ptr1->next;
  }// end of 1st while 

     Binding_Site_Residue_head_ptr=head_atom;

    display_node_list(Binding_Site_Residue_head_ptr);

}//end of 

void Model::SetDifferentAminoAcidTypeInfoForThisModel()
{
  node<Residue> *cursor_aa_ptr;
  node<Aa_Summary> *head_residue_type, *previous_residue_type,*cursor_residue_type;
  previous_residue_type=NULL;
  string aa_list[19]={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU",
                    "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
  int i;

 for(i=0;i<19;i++)
 {
  Aa_Summary *New_aa;
  New_aa=new Aa_Summary;
  char *temp;
  temp=new char;
  strcpy(temp,aa_list[i].c_str());
  New_aa->SetAaName(temp);
  insert_an_node(head_residue_type,previous_residue_type,*New_aa);
 } 

  //node<char*> cursor_ptr;
  cursor_residue_type=head_residue_type;
  while(cursor_residue_type!=NULL)
  {
  cursor_aa_ptr=head_aa_ptr;
  vector<double> temp_v;
  while(cursor_aa_ptr!=NULL)
  {
  if(strncmp(cursor_aa_ptr->data.GetAaName(),cursor_residue_type->data.GetAaName(),3)==0)
  {
   temp_v.push_back(cursor_aa_ptr->data.GetOmega());
  }
  cursor_aa_ptr=cursor_aa_ptr->next;
  }

  cursor_residue_type->data.SetAaNum(temp_v.size());
  cursor_residue_type->data.SetAaAverage(GetModelAaOmegaAverage(temp_v));
  cursor_residue_type->data.SetAaStd(GetModelAaOmegaSd_Of_Ave(temp_v));

  temp_v.clear();
  cursor_residue_type=cursor_residue_type->next;
  }

   aa_sum=head_residue_type;

  //cursor_residue_type=head_residue_type;
  //while(cursor_residue_type!=NULL)
  //{
  // cout<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<p_name<<setw(7)<<model_name<<setw(5)<<model_no
  // <<setw(5)<<cursor_residue_type->data.GetAaName()<<setw(5)
  //    <<cursor_residue_type->data.GetAaNum();
  //cout<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
  //<<setprecision(2)<<setw(7)<<cursor_residue_type->data.GetAaAverage()<<" "<<cursor_residue_type->data.GetAaStd()<<endl;
  //cursor_residue_type=cursor_residue_type->next;
  //} 

}

void Model::CalculateOmegaAngleForAaOfThisModel()
{
  Vector3D V_P_C;
  Vector3D V_Aa_Calpha;
  Vector3D V_Aa_Sc;
  Vector3D V_PC_AaCalpha;
  Vector3D V_AaCalpha_Aa_Sc;

  V_P_C.SetXComponent(GetModelCenter()->data.GetAtom_x());
  V_P_C.SetYComponent(GetModelCenter()->data.GetAtom_y());
  V_P_C.SetZComponent(GetModelCenter()->data.GetAtom_z());


  node<Residue> *cursor_aa_ptr;
  cursor_aa_ptr=head_aa_ptr;

  while(cursor_aa_ptr!=NULL)
  {
    //if(cursor_aa_ptr->data.GetAaCalphaAtom()!=NULL&&cursor_aa_ptr->data.GetAaSideChainGeometricalCenter()!=NULL)
    if(cursor_aa_ptr->data.GetAaCalphaAtom()!=NULL&&cursor_aa_ptr->data.GetAaBetaAtom()!=NULL)
    {
    V_Aa_Calpha.SetXComponent(cursor_aa_ptr->data.GetAaCalphaAtom()->data.GetAtom_x());
    V_Aa_Calpha.SetYComponent(cursor_aa_ptr->data.GetAaCalphaAtom()->data.GetAtom_y());
    V_Aa_Calpha.SetZComponent(cursor_aa_ptr->data.GetAaCalphaAtom()->data.GetAtom_z());

    V_Aa_Sc.SetXComponent(cursor_aa_ptr->data.GetAaBetaAtom()->data.GetAtom_x());
    V_Aa_Sc.SetYComponent(cursor_aa_ptr->data.GetAaBetaAtom()->data.GetAtom_y());
    V_Aa_Sc.SetZComponent(cursor_aa_ptr->data.GetAaBetaAtom()->data.GetAtom_z());

    //V_Aa_Sc.SetXComponent(cursor_aa_ptr->data.GetAaSideChainGeometricalCenter()->data.GetAtom_x());
    //V_Aa_Sc.SetYComponent(cursor_aa_ptr->data.GetAaSideChainGeometricalCenter()->data.GetAtom_y());
    //V_Aa_Sc.SetZComponent(cursor_aa_ptr->data.GetAaSideChainGeometricalCenter()->data.GetAtom_z());
    
    V_PC_AaCalpha=V_Aa_Calpha-V_P_C;      //chain center->residue C alpha atom
    V_AaCalpha_Aa_Sc=V_Aa_Sc-V_Aa_Calpha; //residue C alpha atom->residue side chain geometrical center
    
    cursor_aa_ptr->data.CalculateAndSetOmegaAngleForAa(V_PC_AaCalpha,V_AaCalpha_Aa_Sc);
    }
    cursor_aa_ptr=cursor_aa_ptr->next;
  }

}

void Model::CalculateDistanceBetweenCaOfAaForModel()
{
  node<Residue> *cursor_aa_ptr1;
  node<Residue> *cursor_aa_ptr2;
  node<Residue> *cursor_aa_ptr;

  double distance=0;

  cursor_aa_ptr1=head_aa_ptr;

  while(cursor_aa_ptr1!=NULL)
  {
   cursor_aa_ptr2=head_aa_ptr;
   while(cursor_aa_ptr2!=NULL)
   {
     distance=sqrt(pow((cursor_aa_ptr1->data.GetAaCalphaAtom()->data.GetAtom_x()-
                        cursor_aa_ptr2->data.GetAaCalphaAtom()->data.GetAtom_x()),2)+
                   pow((cursor_aa_ptr1->data.GetAaCalphaAtom()->data.GetAtom_y()-
                        cursor_aa_ptr2->data.GetAaCalphaAtom()->data.GetAtom_y()),2)+
                   pow((cursor_aa_ptr1->data.GetAaCalphaAtom()->data.GetAtom_z()-
                        cursor_aa_ptr2->data.GetAaCalphaAtom()->data.GetAtom_z()),2));
     cursor_aa_ptr2->data.SetAaDistance(distance);
     cursor_aa_ptr1->data.SetAaNeighborAaList(cursor_aa_ptr2->data);
     cursor_aa_ptr2=cursor_aa_ptr2->next;
   }
     cursor_aa_ptr1=cursor_aa_ptr1->next;
  }

  //cursor_aa_ptr1=head_aa_ptr;
  //while(cursor_aa_ptr1!=NULL)
 //{
  //cursor_aa_ptr1->data.SetSortedNeighborAaList();
  //cursor_aa_ptr1=cursor_aa_ptr1->next;
  //}
  

  ofstream protein_aa_neighbor_file(strcat(strtok(p_name,"."),".neighbor"),ios::out);
  cursor_aa_ptr1=head_aa_ptr;
  while(cursor_aa_ptr1!=NULL)
  {
    if(strcmp(cursor_aa_ptr1->data.GetAaName(),"GLY")!=0&&
       cursor_aa_ptr1->data.GetOmega()>=0&&cursor_aa_ptr1->data.GetOmega()<=180)
   { 
   cursor_aa_ptr2=cursor_aa_ptr1->data.GetAaNeighborList();

   node<Residue> *within_5A_residue_list_head,*within_5A_residue_list_previous;
   within_5A_residue_list_head=NULL;
   within_5A_residue_list_previous=NULL;
   int num_of_residue_within_5A;
   num_of_residue_within_5A=0;
 
   while(cursor_aa_ptr2!=NULL)
   {

    if(cursor_aa_ptr2->data.GetAaDistance()>0&&cursor_aa_ptr2->data.GetAaDistance()<=5.0)
     {
      num_of_residue_within_5A=num_of_residue_within_5A+1;
      insert_an_node(within_5A_residue_list_head,within_5A_residue_list_previous,cursor_aa_ptr2->data);
     }
        cursor_aa_ptr2=cursor_aa_ptr2->next;

   }

    protein_aa_neighbor_file<<strtok(p_name,".")<<"\t"<<cursor_aa_ptr1->data.GetAaName()<<"\t"
                            <<num_of_residue_within_5A<<"\t";
    cursor_aa_ptr=within_5A_residue_list_head;
    while(cursor_aa_ptr!=NULL)
   {
    if(cursor_aa_ptr->next!=NULL)
    {
    protein_aa_neighbor_file<<cursor_aa_ptr->data.GetAaName()<<"_"; 
    }
    else
    {
    protein_aa_neighbor_file<<cursor_aa_ptr->data.GetAaName();
    } 
            //<<cursor_aa_ptr->data.GetAaDistance()<<"\t";
        cursor_aa_ptr=cursor_aa_ptr->next;
   }

    protein_aa_neighbor_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)
        <<setw(50-(num_of_residue_within_5A*3+(num_of_residue_within_5A-1)))<<""<<"\t";

    protein_aa_neighbor_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
        <<setprecision(2)<<setw(7)<<cursor_aa_ptr1->data.GetAaAs()
        <<setprecision(2)<<setw(7)<<cursor_aa_ptr1->data.GetAaMs()
        <<setprecision(2)<<setw(7)<<cursor_aa_ptr1->data.GetAaAverageCurvature()
        <<setprecision(2)<<setw(7)<<cursor_aa_ptr1->data.GetOmega()<<endl;
   }

   cursor_aa_ptr1=cursor_aa_ptr1->next;
  }

}

void Model::CalculateDistanceBetweenAtomForModel()
{
  node<Atom> *cursor_atom_ptr1;
  node<Atom> *cursor_atom_ptr2;
  double distance=0;

  cursor_atom_ptr1=head_atom_ptr;

  while(cursor_atom_ptr1!=NULL)
  {
   cursor_atom_ptr2=head_atom_ptr;
   while(cursor_atom_ptr2!=NULL)
   {
     distance=sqrt(pow((cursor_atom_ptr1->data.GetAtom_x()-
                        cursor_atom_ptr2->data.GetAtom_x()),2)+
                   pow((cursor_atom_ptr1->data.GetAtom_y()-
                        cursor_atom_ptr2->data.GetAtom_y()),2)+
                   pow((cursor_atom_ptr1->data.GetAtom_z()-
                        cursor_atom_ptr2->data.GetAtom_z()),2));
     cursor_atom_ptr2->data.SetAtomDistance(distance);
     cursor_atom_ptr1->data.SetAtomNeighborAtomList(cursor_atom_ptr2->data);
     cursor_atom_ptr2=cursor_atom_ptr2->next;
   }
     cursor_atom_ptr1=cursor_atom_ptr1->next;
  }

  cursor_atom_ptr1=head_atom_ptr;
  while(cursor_atom_ptr1!=NULL)
  {
   cursor_atom_ptr1->data.SetSortedNeighborAtomList();
   cursor_atom_ptr1=cursor_atom_ptr1->next;
  }
   
  cursor_atom_ptr1=head_atom_ptr;
 while(cursor_atom_ptr1!=NULL)
  {
   cursor_atom_ptr2=cursor_atom_ptr1->data.GetAtomNeighborList();
 while(cursor_atom_ptr2!=NULL)
 {
    cout<<cursor_atom_ptr1->data.GetAtom_Index()<<" " 
        <<cursor_atom_ptr1->data.GetAtomName()<<" "
        //<<cursor_aa_ptr1->data.GetAaDistance()<<" and "
        <<" and "
        <<cursor_atom_ptr2->data.GetAtom_Index()<<" " 
        <<cursor_atom_ptr2->data.GetAtomName()<<" "
        <<cursor_atom_ptr2->data.GetAtomDistance()<<endl;
 cursor_atom_ptr2=cursor_atom_ptr2->next;
  }
 cursor_atom_ptr1=cursor_atom_ptr1->next;
 }

}

void Model::WriteThisModelSummaryToFile()
{
   node<Residue> *cursor_aa_ptr;
   SetModelConcaveAaList();
   SetModelConvexAaList();
   SetModelSurfaceAaList();
   SetModelBurialAaList();
   SetModelSurfaceConvexAaList();
   SetModelSurfaceConcaveAaList();
 
   node<Aa_Summary> *over_all_aa_head,*over_all_aa_previous;
   over_all_aa_previous=NULL;

    
    ofstream protein_concave_file(strcat(p_name,".concave"),ios::out);
    cursor_aa_ptr=concave_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     protein_concave_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)
     <<cursor_aa_ptr->data.GetAaName()
     <<setw(5)<<cursor_aa_ptr->data.GetAaNum()
     <<setw(5)<<cursor_aa_ptr->data.GetAaChain()
     <<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaMs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAverageCurvature()
     <<setprecision(2)<<setw(9)<<cursor_aa_ptr->data.GetOmega()<<endl;
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while


    ofstream protein_convex_file(strcat(strtok(p_name,"."),".convex"),ios::out);
    cursor_aa_ptr=convex_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     protein_convex_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)
     <<cursor_aa_ptr->data.GetAaName()
     <<setw(5)<<cursor_aa_ptr->data.GetAaNum()
     <<setw(5)<<cursor_aa_ptr->data.GetAaChain()
     <<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaMs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAverageCurvature()
     <<setprecision(2)<<setw(9)<<cursor_aa_ptr->data.GetOmega()<<endl;
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

    ofstream protein_surface_concave_file(strcat(strtok(p_name,"."),".su_ca"),ios::out);
    cursor_aa_ptr=surface_Concave_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     protein_surface_concave_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)
     <<cursor_aa_ptr->data.GetAaName()
     <<setw(5)<<cursor_aa_ptr->data.GetAaNum()
     <<setw(5)<<cursor_aa_ptr->data.GetAaChain()
     <<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaMs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAverageCurvature()
     <<setprecision(2)<<setw(9)<<cursor_aa_ptr->data.GetOmega()<<endl;
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

    ofstream protein_surface_convex_file(strcat(strtok(p_name,"."),".su_cx"),ios::out);
    cursor_aa_ptr=surface_Convex_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     protein_surface_convex_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)
     <<cursor_aa_ptr->data.GetAaName()
     <<setw(5)<<cursor_aa_ptr->data.GetAaNum()
     <<setw(5)<<cursor_aa_ptr->data.GetAaChain()
     <<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaMs()
     <<setprecision(2)<<setw(7)<<cursor_aa_ptr->data.GetAaAverageCurvature()
     <<setprecision(2)<<setw(9)<<cursor_aa_ptr->data.GetOmega()<<endl;
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while
 
    int i;
    char *temp;
    double concave_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double concave_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double convex_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double convex_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double overall_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double overall_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double su_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double su_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double bu_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double bu_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double surface_concave_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double surface_concave_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double surface_convex_av_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double surface_convex_sd_angle[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    temp=new char;

    string aa_list[19]={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU",
                    "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};

    vector<double> convex_value_array;
    vector<double> concave_value_array;
    vector<double> overall_value_array;
    vector<double> su_value_array;
    vector<double> bu_value_array;
    vector<double> surface_convex_value_array;
    vector<double> surface_concave_value_array;

    ofstream protein_file(strcat(strtok(p_name,"."),".summary"),ios::out);
  protein_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"p_n"
    <<setiosflags(ios::right)<<setw(5)<<"aa_t"<<setw(2)<<" ";
  protein_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "
  <<setw(7)<<"ca_av"
  <<setw(7)<<"ca_sd"
  <<setw(7)<<"ca_n"
  <<setw(7)<<"cx_av"
  <<setw(7)<<"cx_sd"
  <<setw(7)<<"cx_n"
  <<setw(7)<<"su_av"
  <<setw(7)<<"su_sd"
  <<setw(7)<<"su_n"
  <<setw(7)<<"bu_av"
  <<setw(7)<<"bu_sd"
  <<setw(7)<<"bu_n"
  <<setw(7)<<"ov_av"
  <<setw(7)<<"ov_sd"
  <<setw(7)<<"ov_n"<<endl;


    ofstream protein_ca_cx_in_surface_file(strcat(strtok(p_name,"."),".ca_cx_in_su"),ios::out);

  /*
  protein_surface_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"p_n"
    <<setiosflags(ios::right)<<setw(5)<<"aa_t"<<setw(2)<<" ";
  protein_surface_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "
  <<setw(7)<<"s_ca_v"
  <<setw(7)<<"s_ca_s"
  <<setw(7)<<"s_ca_n"
  <<setw(7)<<"s_cx_v"
  <<setw(7)<<"s_cx_s"
  <<setw(7)<<"s_cx_n"
  <<setw(7)<<"su_ca_cx"<<endl;
  */

    ofstream protein_surface_burial_file(strcat(strtok(p_name,"."),".bu_su"),ios::out);
  
protein_surface_burial_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<"p_n"
    <<setiosflags(ios::right)<<setw(5)<<"aa_t"<<setw(2)<<" ";

  /*
  protein_surface_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "
  <<setw(7)<<"b_v"
  <<setw(7)<<"b_s"
  <<setw(7)<<"b_n"
  <<setw(7)<<"s__v"
  <<setw(7)<<"s_s"
  <<setw(7)<<"s_n"
  <<setw(7)<<"bu_su"<<endl;
  */


  protein_file<<"p_n\t"<<"aa_t\t"<<setw(2)<<"ca_av\t"<<"ca_sd\t"<<"ca_n\t"<<"cx_av\t"<<"cx_sd\t"
              <<"cx_n\t"<<"su_av\t"<<"su_sd\t"<<"su_n\t"<<"bu_av\t"
              <<"bu_sd\t"<<"bu_n\t"<<"ov_av\t"<<"ov_sd\t"<<"ov_n\t"<<endl;


    Aa_Summary temp_aa;
    double difference_bu_su;
    double difference_ca_cx_in_su;


    for(i=0;i<19;i++)
   {
      
    strcpy(temp,aa_list[i].c_str());
    ofstream aa_file(strcat(temp,".aa"),ios::out|ios::app);
  
    temp_aa.SetAaName(strtok(temp,"."));
    cursor_aa_ptr=concave_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     concave_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      concave_av_angle[i]=GetModelAaOmegaAverage(concave_value_array);
      concave_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(concave_value_array);



    cursor_aa_ptr=convex_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     convex_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      convex_av_angle[i]=GetModelAaOmegaAverage(convex_value_array);
      convex_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(convex_value_array);

    cursor_aa_ptr=surface_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     su_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      su_av_angle[i]=GetModelAaOmegaAverage(su_value_array);
      su_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(su_value_array);


    cursor_aa_ptr=burial_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     bu_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      bu_av_angle[i]=GetModelAaOmegaAverage(bu_value_array);
      bu_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(bu_value_array);

    cursor_aa_ptr=head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     overall_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      overall_av_angle[i]=GetModelAaOmegaAverage(overall_value_array);
      overall_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(overall_value_array);

      temp_aa.SetAaAverage(overall_av_angle[i]);
      temp_aa.SetAaStd(overall_sd_angle[i]);
      temp_aa.SetAaNum(overall_value_array.size());
      insert_an_node(over_all_aa_head,over_all_aa_previous,temp_aa);

    cursor_aa_ptr=surface_Concave_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     surface_concave_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      surface_concave_av_angle[i]=GetModelAaOmegaAverage(surface_concave_value_array);
      surface_concave_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(surface_concave_value_array);

    cursor_aa_ptr=surface_Convex_head_aa_ptr;
    while(cursor_aa_ptr!=NULL)
    {
     if(strncmp(cursor_aa_ptr->data.GetAaName(),temp,3)==0&&cursor_aa_ptr->data.GetOmega()>=0
        &&cursor_aa_ptr->data.GetOmega()<=180)
     {  
     surface_convex_value_array.push_back(cursor_aa_ptr->data.GetOmega());
     }//end if
      cursor_aa_ptr=cursor_aa_ptr->next;
    }//end while

      surface_convex_av_angle[i]=GetModelAaOmegaAverage(surface_convex_value_array);
      surface_convex_sd_angle[i]=GetModelAaOmegaSd_Of_Ave(surface_convex_value_array);

  aa_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<strtok(p_name,".")
    <<setiosflags(ios::right)<<setw(5)<<aa_list[i]<<setw(2)<<" ";
  aa_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "<<setiosflags(ios::fixed)
  <<setprecision(2)<<setw(7)<<concave_av_angle[i]
  <<setprecision(2)<<setw(7)<<concave_sd_angle[i]
  <<setprecision(2)<<setw(7)<<concave_value_array.size()
  <<setprecision(2)<<setw(7)<<convex_av_angle[i]
  <<setprecision(2)<<setw(7)<<convex_sd_angle[i]
  <<setprecision(2)<<setw(7)<<convex_value_array.size()
  <<setprecision(2)<<setw(7)<<su_av_angle[i]
  <<setprecision(2)<<setw(7)<<su_sd_angle[i]
  <<setprecision(2)<<setw(7)<<su_value_array.size()
  <<setprecision(2)<<setw(7)<<bu_av_angle[i]
  <<setprecision(2)<<setw(7)<<bu_sd_angle[i]
  <<setprecision(2)<<setw(7)<<bu_value_array.size()
  <<setprecision(2)<<setw(7)<<overall_av_angle[i]
  <<setprecision(2)<<setw(7)<<overall_sd_angle[i]
  <<setprecision(2)<<setw(7)<<overall_value_array.size()<<endl;

  
  protein_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<strtok(p_name,".")
    <<setiosflags(ios::right)<<setw(5)<<aa_list[i]<<setw(2)<<" ";
  protein_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "<<setiosflags(ios::fixed)
  <<setprecision(2)<<setw(7)<<concave_av_angle[i]
  <<setprecision(2)<<setw(7)<<concave_sd_angle[i]
  <<setprecision(2)<<setw(7)<<concave_value_array.size()
  <<setprecision(2)<<setw(7)<<convex_av_angle[i]
  <<setprecision(2)<<setw(7)<<convex_sd_angle[i]
  <<setprecision(2)<<setw(7)<<convex_value_array.size()
  <<setprecision(2)<<setw(7)<<su_av_angle[i]
  <<setprecision(2)<<setw(7)<<su_sd_angle[i]
  <<setprecision(2)<<setw(7)<<su_value_array.size()
  <<setprecision(2)<<setw(7)<<bu_av_angle[i]
  <<setprecision(2)<<setw(7)<<bu_sd_angle[i]
  <<setprecision(2)<<setw(7)<<bu_value_array.size()
  <<setprecision(2)<<setw(7)<<overall_av_angle[i]
  <<setprecision(2)<<setw(7)<<overall_sd_angle[i]
  <<setprecision(2)<<setw(7)<<overall_value_array.size()<<endl;

 if(surface_concave_value_array.size()>0&&surface_convex_value_array.size()>0&&su_value_array.size()>0)
  { 
   difference_ca_cx_in_su=surface_concave_av_angle[i]-surface_convex_av_angle[i];
  protein_ca_cx_in_surface_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<strtok(p_name,".")
    <<setiosflags(ios::right)<<setw(5)<<aa_list[i]<<setw(2)<<" ";
  protein_ca_cx_in_surface_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "<<setiosflags(ios::fixed)
  <<setprecision(2)<<setw(7)<<surface_concave_av_angle[i]
  <<setprecision(2)<<setw(7)<<surface_concave_sd_angle[i]
  <<setprecision(2)<<setw(7)<<surface_concave_value_array.size()
  <<setprecision(2)<<setw(7)<<surface_convex_av_angle[i]
  <<setprecision(2)<<setw(7)<<surface_convex_sd_angle[i]
  <<setprecision(2)<<setw(7)<<surface_convex_value_array.size()
  <<setprecision(2)<<setw(7)<<difference_ca_cx_in_su<<endl;
 }

 if(su_value_array.size()>0&&bu_value_array.size()>0)
  {
  difference_bu_su=bu_av_angle[i]-su_av_angle[i];
  protein_surface_burial_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(6)<<strtok(p_name,".")
    <<setiosflags(ios::right)<<setw(5)<<aa_list[i]<<setw(2)<<" ";
  protein_surface_burial_file<<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "
  <<setiosflags(ios::fixed)
  <<setprecision(2)<<setw(7)<<bu_av_angle[i]
  <<setprecision(2)<<setw(7)<<bu_sd_angle[i]
  <<setprecision(2)<<setw(7)<<bu_value_array.size()
  <<setprecision(2)<<setw(7)<<su_av_angle[i]
  <<setprecision(2)<<setw(7)<<su_sd_angle[i]
  <<setprecision(2)<<setw(7)<<su_value_array.size()
  <<setprecision(2)<<setw(7)<<difference_bu_su<<endl;
  }

      if(!concave_value_array.empty())
      {
      concave_value_array.clear();
      }

      if(!convex_value_array.empty())
      {
      convex_value_array.clear();
      }

      if(!su_value_array.empty())
      {
      su_value_array.clear();
      }

      if(!bu_value_array.empty())
      {
      bu_value_array.clear();
      }

      if(!overall_value_array.empty())
      {
      overall_value_array.clear();
      }

      if(!surface_concave_value_array.empty())
      {
      surface_concave_value_array.clear();
      }

      if(!surface_convex_value_array.empty())
      {
      surface_convex_value_array.clear();
      }


   }// end for

   node<Aa_Summary> *cursor_aa_summary_ptr1;
   node<Aa_Summary> *cursor_aa_summary_ptr2;
   double difference_for_pair;
   ofstream protein_diff_aa_pair_file(strcat(strtok(p_name,"."),".diff_aa_pair"),ios::out);

   cursor_aa_summary_ptr1=over_all_aa_head;

   while(cursor_aa_summary_ptr1!=NULL)
   {
   
    cursor_aa_summary_ptr2=cursor_aa_summary_ptr1->next;
    while(cursor_aa_summary_ptr2!=NULL)
    {
     if(strcmp(cursor_aa_summary_ptr1->data.GetAaName(),cursor_aa_summary_ptr2->data.GetAaName())!=0&&
               cursor_aa_summary_ptr1->data.GetAaNum()!=0&&cursor_aa_summary_ptr2->data.GetAaNum()!=0)
     {
      difference_for_pair=cursor_aa_summary_ptr1->data.GetAaAverage()-cursor_aa_summary_ptr2->data.GetAaAverage();
      protein_diff_aa_pair_file<<strtok(p_name,".")<<"\t"<<cursor_aa_summary_ptr1->data.GetAaName()<<"_"
                          <<cursor_aa_summary_ptr2->data.GetAaName()<<"\t"
          <<resetiosflags(ios::adjustfield)<<setiosflags(ios::right)<<setw(4)<<" "<<setiosflags(ios::fixed)
          <<setprecision(2)<<setw(7)<<difference_for_pair<<endl;
     }
     cursor_aa_summary_ptr2=cursor_aa_summary_ptr2->next;
    }
     cursor_aa_summary_ptr1=cursor_aa_summary_ptr1->next;
   }

}//end write
 
//void Model::operator =(Model*& protein)
//{
// p_name=new char;
// strcpy(p_name,protein->p_name);
// head_aa_ptr=protein->head_aa_ptr;
//}

void Model::CalculateAndOutputDifferenceInAverageOmegaBetweenBuAndSu()
{


}

void Model::OutputModelSurfaceAalist()
{
   SetModelSurfaceAaList();
   node<Residue> *cursor_aa_ptr;
   cursor_aa_ptr=surface_head_aa_ptr;

   ofstream output(strcat(strtok(p_name,"."),".su_aa"),ios::out);
  
  while(cursor_aa_ptr!=NULL)
  {
  if(strcmp(cursor_aa_ptr->data.GetAaName(),"GLY")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"ASX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"GLX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"PCA")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"UNK")!=0&&
     (cursor_aa_ptr->data.GetOmega()>=0&&cursor_aa_ptr->data.GetOmega()<=180))
   {
   output<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<strtok(p_name,".")
  <<cursor_aa_ptr->data<<endl;
   }
   cursor_aa_ptr=cursor_aa_ptr->next;
  }

}
void Model::OutputModelSurfaceConvexAaList()
{
   SetModelSurfaceConvexAaList();
   node<Residue> *cursor_aa_ptr;
   cursor_aa_ptr=surface_Convex_head_aa_ptr;
   ofstream output(strcat(strtok(p_name,"."),".su_cx_aa"),ios::out);
  
  while(cursor_aa_ptr!=NULL)
  {
  if(strcmp(cursor_aa_ptr->data.GetAaName(),"GLY")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"ASX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"GLX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"PCA")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"UNK")!=0&&
     (cursor_aa_ptr->data.GetOmega()>=0&&cursor_aa_ptr->data.GetOmega()<=180))
   {
   output<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<strtok(p_name,".")
  <<cursor_aa_ptr->data<<endl;
   }
   cursor_aa_ptr=cursor_aa_ptr->next;
  }

}
void Model::OutputModelSurfaceConcaveAaList()
{
   SetModelSurfaceConcaveAaList();
   node<Residue> *cursor_aa_ptr;
   cursor_aa_ptr=surface_Concave_head_aa_ptr;
   ofstream output(strcat(strtok(p_name,"."),".su_ca_aa"),ios::out);
  
  while(cursor_aa_ptr!=NULL)
  {
  if(strcmp(cursor_aa_ptr->data.GetAaName(),"GLY")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"ASX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"GLX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"PCA")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"UNK")!=0&&
     (cursor_aa_ptr->data.GetOmega()>=0&&cursor_aa_ptr->data.GetOmega()<=180))
   {
   output<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<strtok(p_name,".")
  <<cursor_aa_ptr->data<<endl;
   }
   cursor_aa_ptr=cursor_aa_ptr->next;
  }

}
void Model::OutputModelBurialAaList()
{
   SetModelBurialAaList();
   node<Residue> *cursor_aa_ptr;
   cursor_aa_ptr=burial_head_aa_ptr;
   ofstream output(strcat(strtok(p_name,"."),".bu_aa"),ios::out);
  
  while(cursor_aa_ptr!=NULL)
  {
  if(strcmp(cursor_aa_ptr->data.GetAaName(),"GLY")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"ASX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"GLX")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"PCA")!=0&&
     strcmp(cursor_aa_ptr->data.GetAaName(),"UNK")!=0&&
     (cursor_aa_ptr->data.GetOmega()>=0&&cursor_aa_ptr->data.GetOmega()<=180))
   {
   output<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<strtok(p_name,".")
  <<cursor_aa_ptr->data<<endl;
   }
   cursor_aa_ptr=cursor_aa_ptr->next;
  }

}

void Model::OutputModelBackBoneAtom()
{
   node<Atom> *cursor_atom_ptr;
   cursor_atom_ptr=head_atom_ptr;
   
     ofstream output(strcat(strtok(p_name,"."),"_bb.pdb"),ios::out);
   
      while(cursor_atom_ptr!=NULL)
      {
      if(strcmp(cursor_atom_ptr->data.GetAtomName(),"N  ")==0||
         strcmp(cursor_atom_ptr->data.GetAtomName(),"CA ")==0||
         strcmp(cursor_atom_ptr->data.GetAtomName(),"C  ")==0||
         strcmp(cursor_atom_ptr->data.GetAtomName(),"O  ")==0)
      {
       //put<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<strtik(p_name,".")i
      output<<cursor_atom_ptr->data<<endl;
      }

      cursor_atom_ptr=cursor_atom_ptr->next;
      }
   
 }

ostream& operator<<(ostream& os,Model& chain)
{
   /*
   node<Atom> *cursor_atom_ptr;
   cursor_atom_ptr=chain.head_atom_ptr;
   while(cursor_atom_ptr!=NULL)
   {
    os<<cursor_atom_ptr->data<<endl;
   cursor_atom_ptr=cursor_atom_ptr->next;
    }

   node<Residue> *cursor_aa_ptr;
   cursor_aa_ptr=chain.head_aa_ptr;
*/
   node<Residue> *cursor_aa_ptr;
   cursor_aa_ptr=chain.head_aa_ptr;

//os<<chain.model_name<<chain.model_no<<endl;
//<<chain.GetModelCenter()->data<<endl;
 

     
  while(cursor_aa_ptr!=NULL)
  {
  //if(strcmp(cursor_aa_ptr->data.GetAaName(),"GLY")!=0&&
  //   strcmp(cursor_aa_ptr->data.GetAaName(),"ASX")!=0&&
  //  strcmp(cursor_aa_ptr->data.GetAaName(),"GLX")!=0&&
  //   strcmp(cursor_aa_ptr->data.GetAaName(),"PCA")!=0&&
  //   strcmp(cursor_aa_ptr->data.GetAaName(),"UNK")!=0&&
  //   (cursor_aa_ptr->data.GetOmega()>=0&&cursor_aa_ptr->data.GetOmega()<=180))
   //{
   //os<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setw(5)<<chain.p_name
     //<<chain.GetModelCenter()->data<<endl;

   //<<resetiosflags(ios::adjustfield)<<setiosflags(ios::left)<<setiosflags(ios::fixed)
  //<<setprecision(2)<<setw(7)<<chain.GetPrCenter()->data.GetAtom_x()
  //<<setprecision(2)<<setw(7)<<chain.GetPrCenter()->data.GetAtom_y()
  //<<setprecision(2)<<setw(7)<<chain.GetPrCenter()->data.GetAtom_z()
  os<<chain.p_name<<" "<<chain.model_name<<" "<<chain.model_no<<" "<<cursor_aa_ptr->data<<endl;
   //}
   cursor_aa_ptr=cursor_aa_ptr->next;
  }


 return os;

}




