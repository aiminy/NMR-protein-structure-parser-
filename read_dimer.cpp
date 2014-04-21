//@Aimin Yan
//filename:read_dimer.cpp
//to compile: g++ -o read_dimer read_dimer.cpp
//usage: ./readdimer ./*.mmmol
//read dimer protein_rou quartenary file downloaded from ftp://ftp.ebi.ac.uk/pub/databases/msd/pqs/macmol/
///e.g.
//input: ./1a3a_1.mmol
//output:1a3a_1_A.mmol 1a3a_1_C.mmol

#include "cal_tot_asa_of_res.h"
//void geo_cen_of_side_chain3(vector<aa_rou> &w_c3,vector<aa_rou> &w_c2,vector<aa_rou> w_c);
double squ_of_num2(double x);

double squ_of_num2(double x)
{
double y;
y=x*x;
return y;
}

void sum_and_output4(vector<aa_rou> proteins2, vector<aa_rou> proteins1, vector<aa_rou> proteins);

void sum_and_output4(vector<aa_rou> proteins2, vector<aa_rou> proteins1, vector<aa_rou> proteins)
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

  
       //for patch_0_1 
       for(k=0;k<19;k++)
       {  
          angle_total_value[k]=0; 
          num_of_aa[k]=0;
          angle_average_value[k]=0;

          temp1=residue1[k]+".patch_0_1";
          strcpy(temp,temp1.c_str()); 
          cout<<temp1<<endl;
	  ofstream file(temp,ios::out|ios::app);
	  file<<setiosflags(ios::left)<<setw(10)<<"angle"<<"\t"<<"type"<<endl;

	  for(i=0;i<proteins.size();i++)
	  {
	    if(strncmp(proteins[i].res_id,temp,3)==0)
            {
               angle_total_value[k]=angle_total_value[k]+proteins[i].angle;
               angle_total_num=angle_total_num+1;
               num_of_aa[k]=num_of_aa[k]+1;
               file<<setiosflags(ios::left)<<setw(10)<<proteins[i].angle<<"\t"<<proteins[i].res_id<<endl;  
            }//end if
          } //end for

           angle_average_value[k]=angle_total_value[k]/num_of_aa[k];

        file.close();
      } // end of for k=0

        //for patch_1_10  
       for(k=0;k<19;k++)
       {  

          angle_total_value1[k]=0; 
          num_of_aa1[k]=0;
          angle_average_value1[k]=0;
 
          temp1=residue1[k]+".patch_1_10";
          strcpy(temp,temp1.c_str()); 
          cout<<temp1<<endl;
	  ofstream file(temp,ios::out|ios::app);
	  file<<setiosflags(ios::left)<<setw(10)<<"angle"<<"\t"<<"type"<<endl;
	  for(i=0;i<proteins1.size();i++)
	  {
	    if(strncmp(proteins1[i].res_id,temp,3)==0)
            {
               angle_total_value1[k]=angle_total_value1[k]+proteins1[i].angle;
               angle_total_num1=angle_total_num1+1;
               num_of_aa1[k]=num_of_aa1[k]+1;
 
               file<<setiosflags(ios::left)<<setw(10)<<proteins1[i].angle<<"\t"<<proteins1[i].res_id<<endl;  
            }//end if
          }//end for

           angle_average_value1[k]=angle_total_value1[k]/num_of_aa1[k];

        file.close();
      } // end of for k=0

       //for patch_>_10
       for(k=0;k<19;k++)
       {  
          angle_total_value2[k]=0; 
          num_of_aa2[k]=0;
          angle_average_value2[k]=0;
           
          temp1=residue1[k]+".patch_>_10";
          strcpy(temp,temp1.c_str()); 
          cout<<temp1<<endl;
	  ofstream file(temp,ios::out|ios::app);
	  file<<setiosflags(ios::left)<<setw(10)<<"angle"<<"\t"<<"type"<<endl;
	  for(i=0;i<proteins2.size();i++)
	  {
	    if(strncmp(proteins2[i].res_id,temp,3)==0)
            {
               angle_total_value2[k]=angle_total_value2[k]+proteins2[i].angle;
               angle_total_num2=angle_total_num2+1;
               num_of_aa2[k]=num_of_aa2[k]+1;
 
               file<<setiosflags(ios::left)<<setw(10)<<proteins2[i].angle<<"\t"<<proteins2[i].res_id<<endl;  
            }
          }

           angle_average_value2[k]=angle_total_value2[k]/num_of_aa2[k];

        file.close();
      } // end of for k=0

          double pro_aa[19],pro_aa1[19],pro_aa2[19];
          double t;
             
              t=0;

           for(i=0;i<19;i++)
           { 
            pro_aa[i]=0;
            pro_aa1[i]=0;
            pro_aa2[i]=0;
  
            t=num_of_aa[i]/angle_total_num*100;
            pro_aa[i]=t;

            t=num_of_aa1[i]/angle_total_num1*100;
            pro_aa1[i]=t;

            t=num_of_aa2[i]/angle_total_num2*100;
            pro_aa2[i]=t;

           }

          ofstream file1("aa.txt",ios::out|ios::app);
          file1<<"0_1_patch "<<"\t"<<angle_total_num<<"\t"
               <<"1_10_patch"<<"\t"<<angle_total_num1<<"\t"
               <<"10_patch  "<<"\t"<<angle_total_num2<<endl;
          
          
         for(i=0;i<19;i++)
          { 
	  file1<<setiosflags(ios::left)<<setw(10)<<residue1[i]<<setiosflags(ios::fixed)
               <<setiosflags(ios::left)<<setprecision(2)<<setw(5)<<pro_aa[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(5)<<pro_aa1[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(5)<<pro_aa2[i]
               <<endl;
          }

          ofstream file2("omega.txt",ios::out|ios::app);
          file2<<"0_1_patch "<<"\t"<<angle_total_num<<"\t"
               <<"1_10_patch"<<"\t"<<angle_total_num1<<"\t"
               <<"10_patch  "<<"\t"<<angle_total_num2<<endl;
          
          
         for(i=0;i<19;i++)
          { 
	  file2<<setiosflags(ios::left)<<setw(10)<<residue1[i]<<setiosflags(ios::fixed)
               <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<angle_average_value[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<angle_average_value1[i]
               <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<angle_average_value2[i]
               <<endl;
          }

          file1.close(); 

          file2.close(); 
        
}// end of function sum_and_output4

void geo_cen_of_side_chain3(vector<aa_rou> &w_c3,vector <aa_rou> &w_c2,vector<aa_rou> w_c)
{
       
        aa_rou amino;
	amino.res_id       = new char [5];
	
	atom_rou atom_of_res;
	atom_of_res.atom_rou_id = new char [5];
	
	int i,j,k,m;
        //cout<<w_c.size()<<endl;
	
	for(i=0;i<w_c.size();i++)
	{
           //cout<<"test2"<<endl;
	    if(strncmp(w_c[i].res_id,"GLY",3)!=0)
	   {
	    amino.chain_type=w_c[i].chain_type;
	    amino.res_id=w_c[i].res_id;
	    amino.res_num=w_c[i].res_num;
	   
           //for calpha atom
	   k=0;
	   for(j=0;j<w_c[i].atom_rou_str.size();j++)
	   {
             //cout<<"test2"<<endl;
  
            if(strncmp(w_c[i].atom_rou_str[j].atom_rou_id," CA",3)==0)
	    k=k+1;
	   }//end of for(j=0..
	
	   atom_of_res.x=0;
	   atom_of_res.y=0;
	   atom_of_res.z=0;
	   //atom_of_res.asa=0;
	   atom_of_res.b=0;
	   if(!amino.atom_rou_str.empty())
	   {amino.atom_rou_str.clear();}
	  
           //cout<<w_c[0].atom_rou_str[1].atom_rou_id<<endl;
 
	   for(j=0;j<w_c[i].atom_rou_str.size();j++)
	 {
          if(strncmp(w_c[i].atom_rou_str[j].atom_rou_id," CA",3)==0)
	   {
	   atom_of_res.x=atom_of_res.x+w_c[i].atom_rou_str[j].x;
           //cout<<atom_of_res.x<<endl; 
	   atom_of_res.y=atom_of_res.y+w_c[i].atom_rou_str[j].y;
	   atom_of_res.z=atom_of_res.z+w_c[i].atom_rou_str[j].z;
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_rou_str[j].asa;
	   atom_of_res.b=atom_of_res.b+w_c[i].atom_rou_str[j].b;
	   } // end of if
	 } //end of for(j=0..
	
	   atom_of_res.atom_rou_id="CA";
           atom_of_res.x=atom_of_res.x/k;
	   atom_of_res.y=atom_of_res.y/k;
	   atom_of_res.z=atom_of_res.z/k;
	   //atom_of_res.asa=atom_of_res.asa/k;
	   atom_of_res.b=atom_of_res.b/k;
	   amino.atom_rou_str.push_back(atom_of_res);

         // for c_sc
           
	   atom_of_res.x=0;
	   atom_of_res.y=0;
	   atom_of_res.z=0;
	   //atom_of_res.asa=0;
	   atom_of_res.b=0;
	   m=0;
	   for(j=0;j<w_c[i].atom_rou_str.size();j++)
	  {
          if(strncmp(w_c[i].atom_rou_str[j].atom_rou_id," CA",3)!=0&&strncmp(w_c[i].atom_rou_str[j].atom_rou_id," C",3)!=0
	  &&strncmp(w_c[i].atom_rou_str[j].atom_rou_id," N",3)!=0&&strncmp(w_c[i].atom_rou_str[j].atom_rou_id," O",3)!=0)
	   {
	   atom_of_res.x=atom_of_res.x+w_c[i].atom_rou_str[j].x;
	   atom_of_res.y=atom_of_res.y+w_c[i].atom_rou_str[j].y;
	   atom_of_res.z=atom_of_res.z+w_c[i].atom_rou_str[j].z;
	   //atom_of_res.asa=atom_of_res.asa+w_c[i].atom_rou_str[j].asa;
	   atom_of_res.b=atom_of_res.b+w_c[i].atom_rou_str[j].b;
	   m=m+1;
	   }// end of if(stncmp..
	  }// end of for(j=0..
	
           if(m!=0) 
	   {atom_of_res.atom_rou_id="C_SC";
	    atom_of_res.x=atom_of_res.x/m;
	    atom_of_res.y=atom_of_res.y/m;
	    atom_of_res.z=atom_of_res.z/m;
	    //atom_of_res.asa=atom_of_res.asa/m;
	    atom_of_res.b=atom_of_res.b/m;
	    amino.atom_rou_str.push_back(atom_of_res);
           }// end of if
 
           if(k!=0&&m!=0)
           {w_c2.push_back(amino);}//end of if(k!=0&&m!=0) 	 
	 }//end of if("!=GLY)
	
        }// end of for(i<w_c.size()..

         
         double c_x,c_y,c_z;
	 c_x=0;
	 c_y=0;
	 c_z=0;
	 int n;
	 n=0;
	for(i=0;i<w_c2.size();i++)
	{
          for(j=0;j<w_c2[i].atom_rou_str.size();j++)
	 
	{
	c_x=c_x+w_c2[i].atom_rou_str[j].x;
	c_y=c_y+w_c2[i].atom_rou_str[j].y;
	c_z=c_z+w_c2[i].atom_rou_str[j].z;
	n=n+1;
	}// end of for(j=0;
	
	} //end of for(i=0;
	c_x=c_x/n;
	c_y=c_y/n;
	c_z=c_z/n;
 
        
        aa_rou amino_n;
        amino_n.res_id=new char[5];
	
        double a_c_x,a_c_y,a_c_z,t1,t2,t3,ta_x,ta_y,ta_z;
	double sc_a_x,sc_a_y,sc_a_z,angle,ac;
        
	a_c_x=0;a_c_y=0;a_c_z=0;t1=0;t2=0;t3=0;
	sc_a_x=0;sc_a_y=0;sc_a_z=0;
        ta_x=0;ta_y=0;ta_z=0;
	
	for(i=0;i<w_c2.size();i++)
	{
	for(j=0;j<w_c2[i].atom_rou_str.size();j++)
	{
 	 if(strncmp(w_c2[i].atom_rou_str[j].atom_rou_id,"CA",2)==0)
	 {
          t1=sqrt(squ_of_num2(w_c2[i].atom_rou_str[j].x-c_x)+squ_of_num2(w_c2[i].atom_rou_str[j].y-c_y)+
	          squ_of_num2(w_c2[i].atom_rou_str[j].z-c_z));
	  a_c_x=(w_c2[i].atom_rou_str[j].x-c_x)/t1;
	  a_c_y=(w_c2[i].atom_rou_str[j].y-c_y)/t1;
	  a_c_z=(w_c2[i].atom_rou_str[j].z-c_z)/t1;
	  ta_x=w_c2[i].atom_rou_str[j].x;
	  ta_y=w_c2[i].atom_rou_str[j].y;
	  ta_z=w_c2[i].atom_rou_str[j].z;
	  }//end of if  center of structure -> calpha
 	 
	if(strncmp(w_c2[i].atom_rou_str[j].atom_rou_id,"C_SC",4)==0)
	 {
          t2=sqrt(squ_of_num2(w_c2[i].atom_rou_str[j].x-ta_x)+squ_of_num2(w_c2[i].atom_rou_str[j].y-ta_y)+
            squ_of_num2(w_c2[i].atom_rou_str[j].z-ta_z));
	  sc_a_x=(ta_x-w_c2[i].atom_rou_str[j].x)/t2;
	  sc_a_y=(ta_y-w_c2[i].atom_rou_str[j].y)/t2;
	  sc_a_z=(ta_z-w_c2[i].atom_rou_str[j].z)/t2;
	 }//end of if    center of side chain ->calpha
	}// end of for(j=0;j<w_c2[i]
          
        t3=sqrt(squ_of_num2(a_c_x)+squ_of_num2(a_c_y)+squ_of_num2(a_c_z))*sqrt(squ_of_num2(sc_a_x)+squ_of_num2(sc_a_y)+
	        squ_of_num2(sc_a_z));
	ac=acos((a_c_x*sc_a_x+a_c_y*sc_a_y+a_c_z*sc_a_z)/t3);
	angle=180*ac/3.14;

	amino_n.res_num=w_c2[i].res_num;
	amino_n.res_id=w_c2[i].res_id;
	amino_n.chain_type=w_c2[i].chain_type;
	amino_n.angle=angle;

	w_c3.push_back(amino_n);
	}//end of for(i=0;i<w_c2.....)

         //cout<<"this is geo3"<<endl;

}// end of function geo_of_

void read_pdb_file(vector <protein_rou> &protein_rous,char *pdb_file_name)
{

	char dummy[RECORD_MAX];       // temporary C string for reading fields

	ifstream pdb_file2;

	pdb_file2.open(pdb_file_name,ios::in);

        int num_of_protein_rou,i,j;
        num_of_protein_rou=0;
        i=0;  

        atom_rou atom_rou_of_str;
	aa_rou amino;
        protein_rou one_protein_rou;
        //vector <protein_rou> protein_rous;   
    	
	while(pdb_file2>>dummy)
	{
           if(strncmp(dummy,"HEADER",6)==0)
		{
                 //read in a protein_rou   
                 pdb_file2.ignore(56);
                 one_protein_rou.protein_rou_id=new char[12]; 
                 pdb_file2.get(one_protein_rou.protein_rou_id,12);
                 //cout<<one_protein_rou.protein_rou_id<<endl;
                   
                        pdb_file2>>dummy;
                        //cout<<dummy<<endl; 
                 
	 	  while(strncmp(dummy,"ATOM",4)!=0) 
			{
			pdb_file2.getline(dummy,RECORD_MAX);
                        pdb_file2>>dummy;
			}//end of if

		  while(strncmp(dummy,"ATOM",4)==0)
                        {
			pdb_file2.ignore(2);
			pdb_file2.get(dummy,6);
                        //cout<<dummy<<endl;

	                atom_rou_of_str.atom_rou_num=atoi(dummy);
			pdb_file2.ignore(1);
                        
	                atom_rou_of_str.atom_rou_id = new char[6];
			pdb_file2.get(atom_rou_of_str.atom_rou_id,6);
                        //cout<<atom_rou_of_res.atom_rou_id<<endl;   
			
	                amino.res_id= new char [4];
			pdb_file2.get(amino.res_id,4);
                        //cout<<amino.res_id<<endl;	

			pdb_file2.ignore(1);
                        amino.chain_type=new char[5]; 
			pdb_file2.get(amino.chain_type,2);
                        //cout<<amino.chain_type<<endl;

			pdb_file2.get(dummy,5); amino.res_num=atoi(dummy);
			pdb_file2.ignore(4);
                        //cout<<dummy<<endl;

			pdb_file2.get(dummy,9); atom_rou_of_str.x = atof(dummy);

                        //cout<<dummy<<endl;

			pdb_file2.get(dummy,9); atom_rou_of_str.y = atof(dummy);
                        //cout<<dummy<<endl;

			pdb_file2.get(dummy,9); atom_rou_of_str.z = atof(dummy);
                        //cout<<dummy<<endl;

			pdb_file2.get(dummy,7); atom_rou_of_str.occ=atof(dummy);
                        //cout<<dummy<<endl;

			pdb_file2.get(dummy,7); atom_rou_of_str.b = atof(dummy);
                        //cout<<dummy<<endl; 
	           
                        if(one_protein_rou.w_c_rou.empty())
                         {
                        amino.atom_rou_str.push_back(atom_rou_of_str);
                        one_protein_rou.w_c_rou.push_back(amino);
                         }
                         else
                         {
                           if(*amino.res_id==*one_protein_rou.w_c_rou.back().res_id&&
                               amino.res_num==one_protein_rou.w_c_rou.back().res_num&&
                              *amino.chain_type==*one_protein_rou.w_c_rou.back().chain_type)
                             {one_protein_rou.w_c_rou.back().atom_rou_str.push_back(atom_rou_of_str);}
                            else
                             {
                             amino.atom_rou_str.clear(); 
                             amino.atom_rou_str.push_back(atom_rou_of_str);
                             one_protein_rou.w_c_rou.push_back(amino);
                             }  
                          }
                      pdb_file2>>dummy; 
                      //i++;
                      }//end of while ==ATOM
                      
                      //cout<<i<<endl;  
	              protein_rous.push_back(one_protein_rou);
                      one_protein_rou.w_c_rou.clear();
                      amino.atom_rou_str.clear(); 
 
                      //i=0;
		     }//end of if "HEADER"
               else
                 {
			pdb_file2.getline(dummy,RECORD_MAX);
                 } //end of else
       	}  // end of while
	pdb_file2.close();

}

void calculate_roughness_of_residue(vector <protein_rou> &protein_rous)
{

        //cout<<protein_rous.size()<<endl;
        int i,j,k;

        vector<aa_rou> w_c2,w_c2_1,w_c2_2;
        vector<aa_rou> w_c3,w_c3_1,w_c3_2;
        vector<aa_rou> patch,patch1,patch2;

 
        for(i=0;i<protein_rous.size();i++)
        {
          //cout<<protein_rous[i].protein_rou_id<<" has:"<<protein_rous[i].w_c_rou.size()<<" residues "<<endl;
         for(j=0;j<protein_rous[i].w_c_rou.size();j++)
          {
           double roughness_of_residue,temp;
           int num_of_atom_rou_a_residue;
           roughness_of_residue=0;
           num_of_atom_rou_a_residue=0;

          for(k=0;k<protein_rous[i].w_c_rou[j].atom_rou_str.size();k++) 
          {
          temp=0;  
          temp=protein_rous[i].w_c_rou[j].atom_rou_str[k].b-2;
          if(temp>=0)
          {roughness_of_residue=roughness_of_residue+temp;
           num_of_atom_rou_a_residue=num_of_atom_rou_a_residue+1;}
          }

          //cout<<roughness_of_residue<<" "<<num_of_atom_rou_a_residue<<endl;
          if(num_of_atom_rou_a_residue>0)
          {protein_rous[i].w_c_rou[j].roughness=roughness_of_residue/num_of_atom_rou_a_residue;}
          else
          {protein_rous[i].w_c_rou[j].roughness=-1;}

         }  
       }


           ofstream output("results.txt",ios::out);

        for(i=0;i<protein_rous.size();i++)
        {
          //output<<protein_rous[i].protein_rou_id<<" has:"<<protein_rous[i].w_c_rou.size()<<" residues "<<endl;
         for(j=0;j<protein_rous[i].w_c_rou.size();j++)
          {
          if(protein_rous[i].w_c_rou[j].roughness>0&&strcmp(protein_rous[i].w_c_rou[j].res_id,"GLY")!=0)  
          {output<<protein_rous[i].w_c_rou[j].res_id<<" "<<protein_rous[i].w_c_rou[j].roughness<<endl;}
          } 
       }


         //should be wrong, since each patch could not be continuous patch, so we can't calculate 
         //geometrical center for each patch
          
          int t1,t2,t3;
           t1=0;
           t2=0; 
           t3=0;
   
        //extract surface patch based on the roughness of residue
        for(i=0;i<protein_rous.size();i++)
        {
          // cout<<protein_rous[i].protein_rou_id<<" has:"<<protein_rous[i].w_c_rou.size()<<" residues "<<endl;
         
         //cout<<"this protein's patch_0_1 has these residues:"<<endl;
         for(j=0;j<protein_rous[i].w_c_rou.size();j++)
          {
           if(strcmp(protein_rous[i].w_c_rou[j].res_id,"GLY")!=0)
           {
            if(protein_rous[i].w_c_rou[j].roughness>0&&protein_rous[i].w_c_rou[j].roughness<=1)
            {patch.push_back(protein_rous[i].w_c_rou[j]);t1=t1+1;}  //end of patch

            if(protein_rous[i].w_c_rou[j].roughness>1&&protein_rous[i].w_c_rou[j].roughness<=10)
            {patch1.push_back(protein_rous[i].w_c_rou[j]);t2=t2+1;} //end of patch1

            if(protein_rous[i].w_c_rou[j].roughness>10)
            {patch2.push_back(protein_rous[i].w_c_rou[j]);t3=t3+1;} //end of patch2
           }
          }   //end of for(j=0

           //patch 0-1 note: center is changed to 0_1_patch center, However patch could not be continus patch
             geo_cen_of_side_chain3(w_c3,w_c2,patch);

           //patch 1-10 note:center is changed to 1_10_patch center
             geo_cen_of_side_chain3(w_c3_1,w_c2_1,patch1);

           //patch >10  note: center is changed to >_10_patch center
             geo_cen_of_side_chain3(w_c3_2,w_c2_2,patch2);

            
            sum_and_output4(patch2,patch1,patch);
            
           //cout<<patch.size()<<" "<<w_c2.size()<<" "<<w_c3.size()<<endl;
             patch.clear();
             patch1.clear();
             patch2.clear();

             w_c2.clear();
             w_c2_1.clear();
             w_c2_2.clear();

       }//end of for(i=0;

        
	 cout<<w_c3.size()<<endl;
         cout<<w_c2.size()<<endl;
         cout<<"partially"<<t1<<" "<<t2<<" "<<t3<<endl;

        sum_and_output4(w_c3_2,w_c3_1,w_c3);

	return;
 }//end of function of calculate_roughness_of_residue
