#include "cal_tot_asa_of_res.h"
//enum AA {ALA,VAL,PHE,PRO,MET,ILE,LEU,ASP,GLU,LYS,ARG,SER,THR,TYR,HIS,CYS,ASN,GLN,TRP};
//a function to calculate mean of a double array
double calculate_mean_of_double_array(vector<double> w);

double calculate_mean_of_double_array(vector<double> w)
{
   double total,average;
   int i,num;
   
   total=0;average=0;num=0;

   for(i=0;i<w.size();i++)
   {
     total=total+w[i];
     num=num+1;
   }

   average=total/num;

   return average;

}



//a function to calculate standard deviation of a double array 
double calculate_sd_of_double_array(vector<double> w);
double calculate_sd_of_double_array(vector<double> w)
{
 double total,sd,average;
 int i,num;
 
 total=0;sd=0;average=0;num=0;

 average=calculate_mean_of_double_array(w);

 for(i=0;i<w.size();i++)
 {
  total=total+(w[i]-average)*(w[i]-average);
  num=num+1;
  }

 sd=sqrt(total/(num-1));

 return sd;

}



void sum_for_local(vector<aa> &w)
{
  //  AA x;
  //  x=AA(5);
 
   //if(x==ALA)
   //{cout<<"test enum"<<x<<endl;}

    char *temp;
    temp=new char[5];

    string aa_list[19]={"ALA","VAL","PHE","PRO","MET","ILE","LEU","ASP","GLU","LYS",
                    "ARG","SER","THR","TYR","HIS","CYS","ASN","GLN","TRP"};
    double tangle[19],avangle[19],num[19];
    double tangle1[19],avangle1[19],num1[19];
    double tangle2[19],avangle2[19],num2[19];
    double tangle3[19],avangle3[19],num3[19],sd3[19];
    double tangle4[19],avangle4[19],num4[19],sd4[19];
    double d_av_ca[19],d_sd_ca[19];
    double d_av_gcsc[19],d_sd_gcsc[19];
    double avangle5[19],sd5[19];
    double avangle6[19],sd6[19];
    double avangle7[19],sd7[19];
    double avangle8[19],sd8[19];
    double avangle9[19],sd9[19];
    double avangle10[19],sd10[19];

    //cout<<"aa[18] "<<aa[18]<<endl;


    int i,j;

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    tangle[j]=0;
    avangle[j]=0;
    num[j]=0;
   
    for(i=0;i<w.size();i++)
    {
         if(strncmp(w[i].res_id,temp,3)==0)
            { 
              tangle[j]=tangle[j]+w[i].angle;
              num[j]=num[j]+1; 
            }

     }
     avangle[j]=tangle[j]/num[j];
    }

    ofstream file("local.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file<<aa_list[i]<<" "<<avangle[i]<<endl;
   }

    //int i,j;
    //for side chain center to ca * normal vector angle
    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    tangle1[j]=0;
    avangle1[j]=0;
    num1[j]=0;
   
    for(i=0;i<w.size();i++)
    {
         if(strncmp(w[i].res_id,temp,3)==0)
            { 
              tangle1[j]=tangle1[j]+w[i].app_angle;
              num1[j]=num1[j]+1; 
            }

     }
     avangle1[j]=tangle1[j]/num1[j];
    }

    ofstream file1("local_6.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file1<<aa_list[i]<<" "<<avangle1[i]<<endl;
   }

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    tangle2[j]=0;
    avangle2[j]=0;
    num2[j]=0;
   
    for(i=0;i<w.size();i++)
    {
         if(strncmp(w[i].res_id,temp,3)==0)
            { 
              tangle2[j]=tangle2[j]+w[i].angle_atom_3c;
              num2[j]=num2[j]+1; 
            }

     }
     avangle2[j]=tangle2[j]/num2[j];
    }

    ofstream file2("local_atom_3c.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file2<<aa_list[i]<<" "<<avangle2[i]<<endl;
   }

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    tangle3[j]=0;
    avangle3[j]=0;
    num3[j]=0;
   
    for(i=0;i<w.size();i++)
    {
         if(strncmp(w[i].res_id,temp,3)==0&&w[i].angle_atom_3n>=0&&w[i].angle_atom_3n<=180)
            { 
              tangle3[j]=tangle3[j]+w[i].angle_atom_3n;
              num3[j]=num3[j]+1; 
            }

     }
          
      //cout<<aa[j]<<" "<<temp<<" "<<tangle3[j]<<" "<<num3[j]<<endl; 
      avangle3[j]=tangle3[j]/num3[j];
    }

    ofstream file3("local_atom_3n.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file3<<aa_list[i]<<" "<<avangle3[i]<<endl;
   }



    vector<double> temp_array;
    vector<double> d_ca, d_gcsc;

    vector<double> a;
    vector<double> b;
    vector<double> c;
    vector<double> d;

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    ofstream aa_file(temp,ios::out|ios::app);  
    //tangle4[j]=0;
    avangle4[j]=0;
    sd4[j]=0;
    //num4[j]=0;
    d_av_ca[j]=0;
    d_sd_ca[j]=0;
    d_av_gcsc[j]=0;
    d_sd_gcsc[j]=0;
    avangle5[j]=0,sd5[j]=0;
    avangle6[j]=0,sd6[j]=0;
    avangle7[j]=0,sd7[j]=0;
    avangle8[j]=0,sd8[j]=0;

   
    for(i=0;i<w.size();i++)
    {
         if(strncmp(w[i].res_id,temp,3)==0&&w[i].angle_atom_3n>=0&&w[i].angle_atom_3n<=180)
            {
              aa_file<<w[i].res_id<<" "<<w[i].angle_atom_3n<<endl; 
              temp_array.push_back(w[i].angle_atom_3n);
              d_ca.push_back(w[i].distance_calpha_to_3atom_plane);
              d_gcsc.push_back(w[i].distance_scgc_to_3atom_plane);
            }

         if(strncmp(w[i].res_id,temp,3)==0&&w[i].angle_8>=0&&w[i].angle_8<=180)
            { 
              a.push_back(w[i].angle_8);
            }

         if(strncmp(w[i].res_id,temp,3)==0&&w[i].angle_9>=0&&w[i].angle_9<=180)
            { 
              b.push_back(w[i].angle_9);
            }

         if(strncmp(w[i].res_id,temp,3)==0&&w[i].angle_10>=0&&w[i].angle_10<=180)
            { 
              c.push_back(w[i].angle_10);
            }

         if(strncmp(w[i].res_id,temp,3)==0&&w[i].angle_11>=0&&w[i].angle_11<=180)
            { 
              d.push_back(w[i].angle_11);
            }
           

     }

      avangle4[j]=calculate_mean_of_double_array(temp_array);
      sd4[j]=calculate_sd_of_double_array(temp_array);

      avangle5[j]=calculate_mean_of_double_array(a);
      sd5[j]=calculate_sd_of_double_array(a);

      avangle6[j]=calculate_mean_of_double_array(b);
      sd6[j]=calculate_sd_of_double_array(b);

      avangle7[j]=calculate_mean_of_double_array(c);
      sd7[j]=calculate_sd_of_double_array(c);

      avangle8[j]=calculate_mean_of_double_array(d);
      sd8[j]=calculate_sd_of_double_array(d);
      

      d_av_ca[j]=calculate_mean_of_double_array(d_ca);
      d_sd_ca[j]=calculate_sd_of_double_array(d_ca);
       
      d_av_gcsc[j]=calculate_mean_of_double_array(d_gcsc);
      d_sd_gcsc[j]=calculate_sd_of_double_array(d_gcsc);

    temp_array.clear();
    d_ca.clear();
    d_gcsc.clear();
    a.clear();
    b.clear();
    c.clear();
    d.clear();

    }

    ofstream file4("local_atom_3n_sd.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file4<<aa_list[i]<<" "<<avangle4[i]<<" "<<sd4[i]
         <<" "<<avangle5[i]<<" "<<sd5[i]<<" "
         <<" "<<avangle6[i]<<" "<<sd6[i]<<" "
         <<" "<<avangle7[i]<<" "<<sd7[i]<<" "
         <<" "<<avangle8[i]<<" "<<sd8[i]<<" "
         <<" "<<d_av_ca[i]<<" "<<d_sd_ca[i]<<" "
         <<d_av_gcsc[i]<<" "<<d_sd_gcsc[i]<<endl;
   }


    //int counter8;
    //counter8=0;

    vector <aa> w_global_G; // global_angle>90: plane normal point in
    vector <aa> w_global_L; // global_angle<90: plane normal point out

    for(i=0;i<w.size();i++)
    {
       cout<<w[i].res_id<<w[i].global_angle<<endl;
    }

        
    
    for(i=0;i<w.size();i++)
    {
         if(w[i].global_angle>90)
            { 
              w_global_G.push_back(w[i]);
            }
    }

    vector<double> temp_array_1;

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    avangle9[j]=0,sd9[j]=0;
   
    for(i=0;i<w_global_G.size();i++)
    {
         if(strncmp(w_global_G[i].res_id,temp,3)==0&&w_global_G[i].angle_atom_3n>=0
                    &&w_global_G[i].angle_atom_3n<=180)
            { 
              cout<<w_global_G[i].res_id<<w_global_G[i].angle_atom_3n<<endl;     
              temp_array_1.push_back(w_global_G[i].angle_atom_3n);
            }

    }

      avangle9[j]=calculate_mean_of_double_array(temp_array_1);
      sd9[j]=calculate_sd_of_double_array(temp_array_1);
      temp_array_1.clear();
  }


    ofstream file4_G("local_atom_3n_sd_G.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file4_G<<aa_list[i]<<" "<<avangle9[i]<<" "<<sd9[i]<<endl;
   }


    for(i=0;i<w.size();i++)
    {
         if(w[i].global_angle<=90)
            { 
              w_global_L.push_back(w[i]);
            }
    }
   
    vector<double> temp_array_2;

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    avangle10[j]=0,sd10[j]=0;
   
    for(i=0;i<w_global_L.size();i++)
    {
         if(strncmp(w_global_L[i].res_id,temp,3)==0&&w_global_L[i].angle_atom_3n>=0
                    &&w_global_L[i].angle_atom_3n<=180)
            { 
              cout<<w_global_L[i].res_id<<w_global_L[i].angle_atom_3n<<endl;     
              temp_array_2.push_back(w_global_L[i].angle_atom_3n);
            }

    }

      avangle10[j]=calculate_mean_of_double_array(temp_array_2);
      sd10[j]=calculate_sd_of_double_array(temp_array_2);
      temp_array_2.clear();
  }


    ofstream file4_L("local_atom_3n_sd_L.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file4_L<<aa_list[i]<<" "<<avangle10[i]<<" "<<sd10[i]<<endl;
   }

    vector <aa> v1;// sign_calpha>0 && sign_gcsc>0  
    for(i=0;i<w.size();i++)
    {
         if(w[i].sign_calpha>=0&&w[i].sign_gcsc>=0)
            { 
              v1.push_back(w[i]);
            }
    }

    //int counter9;
    //counter9=0;

    vector<aa> v2; 
    for(i=0;i<w.size();i++)
    {
          if(w[i].sign_calpha>=0&&w[i].sign_gcsc<=0)
            { 
              v2.push_back(w[i]);
            }
     }

    //int counter10;
    //counter10=0;
    vector<aa> v3;  
    for(i=0;i<w.size();i++)
    {
         if(w[i].sign_calpha<=0&&w[i].sign_gcsc<=0)
            { 
              v3.push_back(w[i]);
            }
     }

    //int counter11;
    //counter11=0;
 
    vector<aa> v4; 
    for(i=0;i<w.size();i++)
    {
         if(w[i].sign_calpha<=0&&w[i].sign_gcsc>=0)
            { 
              v4.push_back(w[i]);
            }
     }


    cout<<v1.size()<<endl;
    cout<<v2.size()<<endl;
    cout<<v3.size()<<endl;
    cout<<v4.size()<<endl;

    double aa_counter1[19],aa_counter2[19],aa_counter3[19],aa_counter4[19];
    double aa_fre1[19],aa_fre2[19],aa_fre3[19],aa_fre4[19];

    for(j=0;j<19;j++)
   { 
    strcpy(temp,aa_list[j].c_str());
    aa_counter1[j]=0;
    aa_counter2[j]=0;
    aa_counter3[j]=0;
    aa_counter4[j]=0;
    aa_fre1[j]=0;
    aa_fre2[j]=0;
    aa_fre3[j]=0;
    aa_fre4[j]=0;
   
    for(i=0;i<v1.size();i++)
    {
         if(strncmp(v1[i].res_id,temp,3)==0)
           {aa_counter1[j]=aa_counter1[j]+1;}
    }
         
    for(i=0;i<v2.size();i++)
    {
         if(strncmp(v2[i].res_id,temp,3)==0)
           {aa_counter2[j]=aa_counter2[j]+1;}
    }
 
    for(i=0;i<v3.size();i++)
    {
         if(strncmp(v3[i].res_id,temp,3)==0)
           {aa_counter3[j]=aa_counter3[j]+1;}
    }

    for(i=0;i<v4.size();i++)
    {
         if(strncmp(v4[i].res_id,temp,3)==0)
           {aa_counter4[j]=aa_counter4[j]+1;}
    }

    //j aa frequency in v1
    double t_v;
    t_v=v1.size();
    aa_fre1[j]=aa_counter1[j]/t_v;

    // j aa frequency in v2
    t_v=v2.size();
    aa_fre2[j]=aa_counter2[j]/t_v;

    // j aa frequency in v3
    t_v=v3.size();
    aa_fre3[j]=aa_counter3[j]/t_v;

    // j aa frequency in v4
    t_v=v4.size();
    aa_fre4[j]=aa_counter4[j]/t_v;
    
   }//end for(j=0;j<19;j++)


    ofstream file5("four_case_fre.txt",ios::out|ios::app); 

   for(i=0;i<19;i++)
   {
    file5<<aa_list[i]<<" "<<aa_fre1[i]<<" "<<aa_fre2[i]<<" "<<aa_fre3[i]<<" "<<aa_fre4[i]<<endl;
   }


  
}// end of function sum_for_local
