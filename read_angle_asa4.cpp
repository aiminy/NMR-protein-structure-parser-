#define PDB_INDEX  1
#define RECORD_MAX  90
#define FIELD_MAX 12
#include <ctype.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>

using namespace std;

typedef struct {
	int    num_of_res;	// total number of residues
	char   **res_id;	// array (#atom x 1) residue type
	char   *chain_id;	// array (#atom x 1) chain identifier
	int    *res_num;	// array (#atom x 1) residue number
	double *asa;		// array (#atom x 1) residue asa
} asa_struct;

typedef struct {
        int    num_of_res;      // total number of residues
        char   **res_id;        // array (#atom x 1) residue type
        char   *chain_id;       // array (#atom x 1) chain identifier
        int    *res_num;        // array (#atom x 1) residue number
	double *asa;
        double *delt_asa;       // array (#atom x 1) residue asa
} temp_struct;

typedef struct {
        int    num_of_res;       //num_of_residue for this protein
        char   **res_id;         // array (#atom x 1) residue type
        char   *chain_id;        // array (#atom x 1) chain identifier
        int    *res_num;         // array (#atom x 1) residue number
        double *rj;              // array (#atom x 1) residue asa
} aa_struct;

typedef struct {
        int    num_of_res;        // total number of residues
        char   **res_id;          // array (#atom x 1) residue type
        char   *chain_id;         // array (#atom x 1) chain identifier
        int    *res_num;          // array (#atom x 1) residue number
        double *angle;            // array (#atom x 1) residue angle
	double *asa;
	double *delta_asa;        //array(#residuex1) residue delta_asa
} angle_asa_struct;

void read_asa_file(asa_struct &mol, char *asa_file_name);
void read_aa_file(aa_struct &aa, char *aa_file_name);

main (int argc, char *argv[]) {

	asa_struct  mol;
	aa_struct aa;
	temp_struct tmp;
	angle_asa_struct rj_asa;
        int    i,j,k,l;
        char out_name1[13],out_name2[13];
	char p_n[8];
	vector <char> chain;
	vector <asa_struct> mol_a;
	vector <aa_struct>  angle_a,angle_w;
	char *c_a;
	c_a=new char[1];
	
	if (argc < PDB_INDEX+1) 
	{
 	 cout <<"Usage: prompt> "<<argv[0]<<" file.rsa\n";
	 exit(1);
	}
	
	read_asa_file(mol, argv[PDB_INDEX]);

        tmp.num_of_res = mol.num_of_res; 
        tmp.res_id     = new char *[tmp.num_of_res];

   	for (i=0; i<tmp.num_of_res; ++i) 
	{
           tmp.res_id[i]= new char [5];
   	}
   	
	tmp.chain_id = new char  [tmp.num_of_res];
   	tmp.res_num  = new int    [tmp.num_of_res];
   	tmp.delt_asa = new double [tmp.num_of_res];

	tmp.asa=new double[tmp.num_of_res];
	
        cout<<mol.chain_id[0]<<endl;

	j=0;
	chain.push_back(mol.chain_id[0]);
	for(i=1;i<mol.num_of_res;i++)
	{
         if(mol.chain_id[i]!=chain[j])
         {chain.push_back(mol.chain_id[i]);j++;}
        }
        
	mol_a.push_back(mol);
	for(j=0;j<chain.size();j++)
        {
	cout<<chain[j]<<endl;
	strcpy(p_n,strtok(argv[1],"."));
        cout<<p_n<<endl;  
	c_a[0]=chain[j];
        strcat(p_n,"_");
	strcat(p_n,c_a);
        cout<<p_n<<endl; 
	strcat(p_n,".rsa");
	ifstream file(p_n,ios::in);
	read_asa_file(mol,p_n);
	mol_a.push_back(mol);
	file.close();
	} //end for

//cout<<mol_a[0].num_of_res<<endl;
        for(j=0;j<chain.size();j++)
       {
	strcpy(p_n,strtok(argv[1],"."));
	c_a[0]=chain[j];
	strcat(p_n,c_a);
	strcat(p_n,"_angle.txt");
	ifstream file(p_n,ios::in);
	read_aa_file(aa,p_n);
	angle_a.push_back(aa);
	file.close();
       } //end for

//cout<<angle_a[0].num_of_res<<angle_a[1].num_of_res<<endl;

         k=0;

	for(i=1;i<mol_a.size();i++)
	{
  	 for(j=0;j<mol_a[i].num_of_res;j++)
   	 {
     	  if((strcmp(mol_a[i].res_id[j],mol_a[0].res_id[k])==0)
             &&mol_a[i].chain_id[j]==mol_a[0].chain_id[k]
             &&mol_a[i].res_num[j]==mol_a[0].res_num[k])
	   {
	    strcpy(tmp.res_id[k],mol_a[0].res_id[k]);
	    tmp.chain_id[k]=mol_a[0].chain_id[k];
	    tmp.res_num[k]=mol_a[0].res_num[k];
	    tmp.asa[k]=mol_a[0].asa[k];
	    tmp.delt_asa[k]=mol_a[i].asa[j]-mol_a[0].asa[k];
	   }
	  k++;
         }//end of for
       }
//cout<<tmp.res_id[58]<<tmp.chain_id[58]<<tmp.res_num[58]<<tmp.delt_asa[58]
//<<endl;
   
        rj_asa.num_of_res=tmp.num_of_res;
        rj_asa.res_id    = new char *[rj_asa.num_of_res];

	for (i=0; i<rj_asa.num_of_res; ++i) 
	{
        rj_asa.res_id[i]    = new char [5];
	}

	rj_asa.chain_id     = new char  [rj_asa.num_of_res];

	rj_asa.res_num = new int [rj_asa.num_of_res];
	rj_asa.angle=new double [rj_asa.num_of_res];
	rj_asa.delta_asa= new double [rj_asa.num_of_res];
	rj_asa.asa=new double[rj_asa.num_of_res];

	l=0;
	for(i=0;i<angle_a.size();i++)
	{
  	for(j=0;j<angle_a[i].num_of_res;j++)
   	{
	for(k=0;k<tmp.num_of_res;k++)
	{
     	if((strcmp(angle_a[i].res_id[j],tmp.res_id[k])==0)
        &&angle_a[i].chain_id[j]==tmp.chain_id[k]
        &&angle_a[i].res_num[j]==tmp.res_num[k])
        {
        strcpy(rj_asa.res_id[l],tmp.res_id[k]);
        rj_asa.chain_id[l]=tmp.chain_id[k];
        rj_asa.res_num[l]=tmp.res_num[k];
	rj_asa.angle[l]=angle_a[i].rj[j];
        rj_asa.delta_asa[l]=tmp.delt_asa[k];
	rj_asa.asa[l]=tmp.asa[k];
	l++;
        }//end of if
	}//end of for k
   	}//end of for for j
	}//end of for i


	rj_asa.num_of_res=l;
//	cout<<tmp.num_of_res<<k<<endl;
//	cout<<rj_asa.num_of_res<<endl;
	strcpy(p_n,strtok(argv[1],"."));
        cout<<"get "<<p_n<<endl;
        strcpy(out_name1,p_n);
	strcat(out_name1,".inter1");
	// strcat(out_name1,".txt");
	strcpy(out_name2,p_n);
	strcat(out_name2,".noninter1");
	// strcat(out_name2,".txt");
 // cout<<out_name1<<out_name2<<endl;
  	ofstream aa_inter,aa_noninter;
  	aa_inter.open(out_name1,ios::out);
  	aa_noninter.open(out_name2,ios::out);

	for(i=0;i<rj_asa.num_of_res;i++)
	{
	if(rj_asa.delta_asa[i]>1)
	{
	aa_inter<<setiosflags(ios::left)<<setw(4)<<rj_asa.res_id[i]
        <<setw(1)<<rj_asa.chain_id[i]
        <<setiosflags(ios::right)<<setw(5)<<rj_asa.res_num[i]
        <<" "<<setiosflags(ios::fixed)
        <<setprecision(3)<<setw(8)<<rj_asa.angle[i]
        <<" "<<setiosflags(ios::right)
	<<setprecision(3)<<setw(8)<<rj_asa.asa[i]<<" " 
        <<setiosflags(ios::right)
        <<setprecision(3)<<setw(8)<<rj_asa.delta_asa[i]<<endl;
	}
	else
	{
	aa_noninter<<setiosflags(ios::left)<<setw(4)<<rj_asa.res_id[i]
        <<setw(1)<<rj_asa.chain_id[i]
        <<setiosflags(ios::right)<<setw(5)<<rj_asa.res_num[i]
        <<" "<<setiosflags(ios::fixed)
        <<setprecision(3)<<setw(8)<<rj_asa.angle[i]
        <<" "<<setiosflags(ios::right)
	<<setprecision(3)<<setw(8)<<rj_asa.asa[i]<<" " 
        <<setiosflags(ios::right)
        <<setprecision(3)<<setw(8)<<rj_asa.delta_asa[i]<<endl;
	}

	}
	aa_inter.close();
	aa_noninter.close(); 

	return(0);
}

void read_asa_file(asa_struct &mol, char *asa_file_name)
{
	ifstream         asa_file,asa_file2;          // stream name for input asa file
	char    dummy[RECORD_MAX];          // temporary C string for reading fields
	register int  i;		    // fast access iterator

	mol.num_of_res=0;
	
	asa_file.open(asa_file_name,std::ios::in);

	if (asa_file.bad()) 
	{
	cerr <<"Error: Unable to open "<<asa_file_name<<endl;
	exit (2);
	}
	while ((asa_file >> dummy)) 
	{
               // cout<<dummy<<endl;
		if (strncmp(dummy,"RES",3)==0) 
		{
		++mol.num_of_res;
                //cout<<mol.num_of_res<<endl;
		}
		asa_file.getline(dummy,RECORD_MAX);
	}
	asa_file.close();
        
       // cout<<"ok"<<mol.num_of_res<<endl;
	if (mol.num_of_res==0) 
	{
	cerr <<mol.num_of_res<<"No usable residues in this asa file.\n";
	exit (1);
	}

	mol.res_id       = new char *[mol.num_of_res];
	for (i=0; i<mol.num_of_res; ++i) 
	{
	mol.res_id[i] = new char [5];
	}
	mol.chain_id = new char  [mol.num_of_res];

	mol.res_num = new int    [mol.num_of_res];
	mol.asa       = new double [mol.num_of_res];
	i=0;

	asa_file2.open(asa_file_name,std::ios::in);
        //cout<<"ok is before while"<<endl;
	while (asa_file2>>dummy){
               // cout<<"ok is after while"<<endl;
		if (strncmp(dummy,"RES",3)!=0) 
		{
	        //cout<<"not res"<<endl;
		asa_file2.getline(dummy,RECORD_MAX);
		}
		else {
			asa_file2.ignore(1);
			asa_file2.get(mol.res_id[i],4);
			asa_file2.ignore(1);
			asa_file2.get(mol.chain_id[i]);

                        //cout<<"chain"<<mol.chain_id[i]<<endl;

			if(isspace(asa_file2.peek()))
			{ 
			asa_file2.ignore(1);
			asa_file2.get(dummy,4); mol.res_num[i]=atoi(dummy);
			}
			else
			{
			asa_file2.get(dummy,5); mol.res_num[i]=atoi(dummy);
			}
			asa_file2.ignore(23);
			asa_file2.get(dummy,6); mol.asa[i] =atof(dummy);
			i++;
			asa_file2.getline(dummy,RECORD_MAX);
		     }
	}
	asa_file2.close();	

	return;
}

void read_aa_file(aa_struct &aa, char *aa_file_name)
{
	std::ifstream      aa_file;     // stream name for input aa file
	char    dummy[RECORD_MAX];       // temporary C string for reading fields
	register int  i;		 // fast access iterator

	aa.num_of_res=0;
	
	aa_file.open(aa_file_name,std::ios::in);
	if (aa_file.bad()) {
		std::cerr <<"Error: Unable to open "<<aa_file_name<<std::endl;
		exit (2);
	}
	while (aa_file.peek()!=-1) {
		aa.num_of_res++;
		aa_file.getline(dummy,RECORD_MAX);
	}
	aa_file.close();
	if (aa.num_of_res==0) {
		std::cerr << "No usable residues in this aa file.\n";
		exit (1);
	}

	aa.res_id       = new char *[aa.num_of_res];
	for (i=0; i<aa.num_of_res; ++i) {
		aa.res_id[i]    = new char [5];
	}
	aa.chain_id     = new char  [aa.num_of_res];

	aa.res_num = new int    [aa.num_of_res];
	aa.rj       = new double [aa.num_of_res];

	i=0;
	aa_file.open(aa_file_name,std::ios::in);
      	while (aa_file.peek()!=-1) {
       	 aa_file.get(aa.res_id[i],4);
         aa_file.get();
         aa_file.get(aa.chain_id[i]);
         aa_file.get();
         aa_file>>dummy;aa.res_num[i]=std::atoi(dummy);
         aa_file.get();
         aa_file.get(dummy,8); aa.rj[i] = std::atof(dummy);
         i++;
         aa_file.getline(dummy,RECORD_MAX);
         }
	aa_file.close();	

	return;
}
