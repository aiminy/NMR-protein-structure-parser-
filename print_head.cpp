
#include "cal_tot_asa_of_res.h"

void print_head(ofstream &num_file, ofstream &angle_file)
{

        //ofstream file("num_file.txt",ios::out);
        num_file<<setiosflags(ios::left)<<setw(12)<<"p_name"
            <<setiosflags(ios::left)<<setw(8)<<"aa_name"
            <<setiosflags(ios::left)<<setw(6)<<"toal#"
            <<setiosflags(ios::left)<<setw(6)<<"#_10"
            <<setiosflags(ios::left)<<setw(6)<<"#_20"
            <<setiosflags(ios::left)<<setw(6)<<"#_30"
            <<setiosflags(ios::left)<<setw(6)<<"#_40"
            <<setiosflags(ios::left)<<setw(6)<<"#_50"
            <<setiosflags(ios::left)<<setw(6)<<"#_60"
            <<setiosflags(ios::left)<<setw(6)<<"#_70"
            <<setiosflags(ios::left)<<setw(6)<<"#_80"
            <<setiosflags(ios::left)<<setw(6)<<"#_90"
            <<setiosflags(ios::left)<<setw(6)<<"#_100"<<endl;
        //file.close();

        //ofstream file1("angle_file.txt",ios::out|ios::app);
        angle_file<<setiosflags(ios::left)<<setw(12)<<"p_name"
             <<setiosflags(ios::left)<<setw(8)<<"aa_name"
             <<setiosflags(ios::left)<<setw(6)<<"toalA"
             <<setiosflags(ios::left)<<setw(6)<<"A_10"
             <<setiosflags(ios::left)<<setw(6)<<"A_20"
             <<setiosflags(ios::left)<<setw(6)<<"A_30"
             <<setiosflags(ios::left)<<setw(6)<<"A_40"
             <<setiosflags(ios::left)<<setw(6)<<"A_50"
             <<setiosflags(ios::left)<<setw(6)<<"A_60"
             <<setiosflags(ios::left)<<setw(6)<<"A_70"
             <<setiosflags(ios::left)<<setw(6)<<"A_80"
             <<setiosflags(ios::left)<<setw(6)<<"A_90"
             <<setiosflags(ios::left)<<setw(6)<<"A_100"<<endl;
        //file1.close();

}
