//#define global_angle 0;
#ifndef CAL_TOT_ASA_OF_RES_H 
#define CAL_TOT_ASA_OF_RES_H

#define PDB_INDEX  1
#define RECORD_MAX  90
#define FIELD_MAX 12

#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <ctype.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <math.h>
//#include "point.h"
#include "linked_list.h"
#include "aa.h"
#include "protein.h"


using namespace std;

typedef struct{
int atom_rou_num;
char *atom_rou_id; //array of character
double x;
double y;
double z;
double occ;
double b;
} atom_rou;

typedef struct{

double asa;
double angle;
double roughness;
int res_num;
char *chain_type;
char *res_id;
vector <atom_rou> atom_rou_str;
} aa_rou;


typedef struct {

char *protein_rou_id;
vector <aa_rou> w_c_rou;
} protein_rou;


typedef struct {
int res_num;
char chain_type;
char *res_id;
double angle;
double asa;
double ave_x;
double ave_y;
double ave_z;
double ave_b;
double delta_asa;
double distance;
double roughness;
} residue;

typedef struct {
char *protein_name;
vector <residue> residue_str;
} protein;

typedef struct {
int res_num;
char chain_type;
char *res_id;
double angle;
double asa;
double angle_t;
double angle_chain_c;
double angle_chain_t;
double asa_chain;
double delta_asa;
double distance;

} naa;

typedef struct {

	int    num_of_atoms;	     // total number of atoms
	int    num_of_res;	     // total number of residues
	int    *atom_num;	     // array (#atom x 1) atom number
	char   **atom_id;	     // array (#atom x 1) atom label
	char   **res_id;	     // array (#atom x 1) residue type
	char   *chain_id;	     // array (#atom x 1) chain identifier
	int    *res_num;	     // array (#atom x 1) residue number
	double *x;		     // array (#atom x 1) x coordinate
	double *y;		     // array (#atom x 1) y coordinate
	double *z;		     // array (#atom x 1) z coordinate
	double *asa;		     // array (#atom x 1) atomic asa

} coord_struct;

typedef struct{
int atomID;
char *atom_id; //array of character
double x;
double y;
double z;
double asa;
double b;
int res_num;
char chain_type;
char *res_id;
} atom;

typedef struct{
int ID1;
int ID2;
double d;
} atom_pair;

typedef struct {
int res_num;
char chain_type;
char *res_id;
vector <atom> atom_str; //array of atom
 
double angle;        // two nearest atoms  
double app_angle;   //  6 nearest atoms
double angle_atom_3c; // triangle from three closest atom: center->calpha . sc->calpha 
double angle_atom_3n;  // normal . sc->calpha
double distance_calpha_to_3atom_plane;// ca distance to plane 
double distance_scgc_to_3atom_plane;  // scgc distance to plane
double angle_8;
double angle_9;
double angle_10;
double angle_11;
double sign_calpha;
double sign_gcsc;
double global_angle;
//vector <Point> point_array; 
} aa;

typedef struct {
int res_num;
char chain_type;
char *res_id;
double tot_atom_asa;
} naa2;

typedef struct {
char *res_id;
double asa_1;
double asa_2;
double asa_3;
double asa_4;
double asa_5;
double asa_6;
double asa_7;
double asa_8;
double asa_9;
double asa_10;
double asa_11;
double asa_12;
double asa_13;
double asa_14;
}aa_asa;

typedef struct{
char *res_id;
char *r_atom;
} std_aa;

double squ_of_num(double x);
void ter_atom_of_side_chain(vector<naa> &w_c3,vector <aa> &w_c2,vector<aa> w_c,char asa[],char *flag);
void ter_atom_of_side_chain2(vector<naa> &w_c3,vector <aa> &w_c2,vector<aa> w_c);
void produce_std_aa(vector<std_aa>&st_aa); 
void geo_cen_of_side_chain(vector<naa> &w_c3,vector<aa> &w_c2,vector<aa> w_c,char asa[],char *flag);
void geo_cen_of_side_chain2(vector<naa> &w_c3,vector<aa> &w_c2,vector<aa> w_c);
void cal_tot_asa_of_res(vector<naa2>&w_c4,vector<aa>w_c);
void read_standard_data(vector<aa_asa> &std_asa);
void read_asa_file(vector <aa> &w_c,coord_struct &asa, char *asa_file_name);
void cal_n_angle_0_180(char asa[],vector<naa>w_c3,char *flag);
void print_head(ofstream &num_file,ofstream &angle_file);
void cal_n_angle_0_180_2(char p_n[],vector<naa>w_c3,double a,double b, ofstream &file1, ofstream &file2);
void sum_and_output(vector<protein> protein_inter);
void sum_count_output(vector<protein> protein_inter);
void sum_and_output3(vector<protein> protein_inter);
void geo_cen_of_res(vector<naa> &w_c3,vector <aa> &w_c2,vector<residue> &a_w_c,vector<aa> w_c);
void res_alpha_atom(vector<naa> &w_c3,vector <aa> &w_c2,vector<residue> &a_w_c,vector<aa> w_c);
void cal_radius_of_gyration(vector<naa> &w_c3,vector <aa> &w_c2,vector<residue> &a_w_c,vector<aa> w_c);
void output_gene_info_a_protein(vector<protein> whole_protein);
void sum_based_on_radius(vector<protein> whole_protein);
void inter_atom_distance_matrix(vector<aa> &w_c,vector<aa> &sp);
void read_pdb_file(vector <protein_rou> &protein_rous, char *pdb_file_name);
void calculate_roughness_of_residue(vector <protein_rou> &protein_rous);
void sum_and_output5(vector<protein> whole_protein);
void sum_for_local(vector<aa> &w);
void calculate_three_nearest_neighbors(vector<aa> &w_c,vector<aa> &sp);
void regenerate_pdb_file(vector<aa> w);
void output_global_angle(vector<aa> w);

//template <class Item>
void read_based_atom(Protein*& protein,char *file);
void inter_atom_distance_matrix2(node<Residue> *&head_aa);

//#include "read_based_atom.cpp"
#endif
