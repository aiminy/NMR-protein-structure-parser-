#include "cal_tot_asa_of_res.h"
#include "tmatrix.h"
#include "Vector3D.h"

template <class T>

int partition(vector<double> &, int, int);
void quicksort(vector<double> &, int, int);
//void assign_atom_pair(atom_pair a, atom_pair b);


//void assign_atom_pair(atom_pair a,atom_pair b)
//{
//     a.ID1=b.ID1;
//     a.ID2=b.ID2;
//     a.d=b.d;
//}

//Function to determine the partitions
// partitions the array and returns the middle index (subscript)
int partition(vector <atom_pair> &array, int top, int bottom)
{
     atom_pair x1;
     //assign_atom_pair(x1,array[top]);

     x1=array[top];    
     int i = top - 1;
     int j = bottom + 1;
     atom_pair temp;

     do
     {
           do      
           {
           j--;
           }while (x1.d>array[j].d);

          do  
          {
          i++;
          } while (x1.d <array[i].d);

          if (i < j)
         { 
          temp = array[i];    // switch elements at positions i and j
          //assign_atom_pair(temp,array[i]); 
          array[i] = array[j];
          //assign_atom_pair(array[i],array[j]); 
          array[j] = temp;
          //assign_atom_pair(array[j],temp); 
         }

     }while (i < j);
     
     return j;           // returns middle index 
}

void quicksort(vector <atom_pair> &array, int top, int bottom)
{
      // top = subscript of beginning of vector being considered
      // bottom = subscript of end of vector being considered
      // this process uses recursion - the process of calling itself
     int middle;
     if (top < bottom)
    {
          middle = partition(array, top, bottom);
          quicksort(array, top, middle);   // sort top partition
          quicksort(array, middle+1, bottom);    // sort bottom partition
     }
     return;
}


void inter_atom_distance_matrix(vector<aa> &w_c,vector<aa> &sp)
{
        int i,j,rows,cols,num_total_atom,num_of_calpha;
            

         num_total_atom=0;
	
	for(i=0;i<w_c.size();i++)
	 {
          num_total_atom=num_total_atom+w_c[i].atom_str.size();
         }

         rows=cols=num_total_atom;
  
       
        vector <atom> chain_all_atom;
 
        for(i=0;i<w_c.size();i++)
         {
          for(j=0;j<w_c[i].atom_str.size();j++)
           {
            w_c[i].atom_str[j].res_num=w_c[i].res_num;
            w_c[i].atom_str[j].chain_type=w_c[i].chain_type;
            w_c[i].atom_str[j].res_id=new char[5];
            w_c[i].atom_str[j].res_id=w_c[i].res_id;
            chain_all_atom.push_back(w_c[i].atom_str[j]);
           }
         }

         vector<atom> surface_calpha_atom;
         int k=1;
         for(i=0;i<chain_all_atom.size();i++)
          {
           if((strcmp(chain_all_atom[i].atom_id,"CA ")==0)&&chain_all_atom[i].asa>5&&
               strncmp(chain_all_atom[i].res_id,"GLY",3)!=0)
             {
              chain_all_atom[i].atomID=k; 
              surface_calpha_atom.push_back(chain_all_atom[i]);
              k++;
             }//end if
          }//end for
 
           /*
           for(i=0;i<surface_calpha_atom.size();i++)
           {
              cout<<setiosflags(ios::left)<<setw(5)<<surface_calpha_atom[i].atomID
                  <<setiosflags(ios::left)<<setw(5)<<surface_calpha_atom[i].res_num
              <<setiosflags(ios::left)<<setw(5)<<surface_calpha_atom[i].chain_type
              <<setiosflags(ios::left)<<setw(5)<<surface_calpha_atom[i].res_id
              <<setiosflags(ios::left)<<setw(5)<<surface_calpha_atom[i].atom_id
              <<setiosflags(ios::fixed)
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<surface_calpha_atom[i].x
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<surface_calpha_atom[i].y
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<surface_calpha_atom[i].z
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<surface_calpha_atom[i].asa<<endl;
          }
           */

          /* 
           for(i=0;i<chain_all_atom.size();i++)
           {
              cout<<setiosflags(ios::left)<<setw(5)<<chain_all_atom[i].atomID
                  <<setiosflags(ios::left)<<setw(5)<<chain_all_atom[i].res_num
              <<setiosflags(ios::left)<<setw(5)<<chain_all_atom[i].chain_type
              <<setiosflags(ios::left)<<setw(5)<<chain_all_atom[i].res_id
              <<setiosflags(ios::left)<<setw(5)<<chain_all_atom[i].atom_id
              <<setiosflags(ios::fixed)
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<chain_all_atom[i].x
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<chain_all_atom[i].y
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<chain_all_atom[i].z
              <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<chain_all_atom[i].asa<<endl;
          }
           */

          
          
         rows=cols=surface_calpha_atom.size();

        
         tmatrix<atom_pair> distance_matrix(rows,cols);
     
         //atom_pair one_atom_pair;
         double t;
         
         t=0; 
         for(i=0;i<surface_calpha_atom.size();i++)
         {
           for(j=0;j<surface_calpha_atom.size();j++)
            {
             distance_matrix[i][j].ID1=surface_calpha_atom[i].atomID;
             distance_matrix[i][j].ID2=surface_calpha_atom[j].atomID;

             distance_matrix[i][j].d=sqrt(pow((surface_calpha_atom[i].x-surface_calpha_atom[j].x),2)+
                    pow((surface_calpha_atom[i].y-surface_calpha_atom[j].y),2)+
                    pow((surface_calpha_atom[i].z-surface_calpha_atom[j].z),2));
            }
         }

        //cout<<distance_matrix[0][1].distance<<endl;

        
         vector<atom_pair> min;
         
         //double t_min;

         //int k; 
         
        for(i=0;i<rows;i++)
        {
         for(j=0;j<cols;j++)
          {
           if(distance_matrix[i][j].d!=0)
           {
           min.push_back(distance_matrix[i][j]);
          }
          }
       } 

        //vector<double> temp;  
        //for(i=0;i<4;i++)
        //{ 
        // temp.push_back(min[i]);
        //}

         
        //qsort(temp);

        /* 
        cout<<"before sorting"<<endl;
        for(i=0;i<min.size();i++)
        { 
         cout<<"("<<min[i].ID1<<","<<min[i].ID2<<","<<min[i].d<<")"<<endl;
        }
        */

 
        quicksort(min, 0, min.size());
        
        /*  
        cout<<"After sorting"<<endl;
        for(i=0;i<min.size();i++)
        {
         cout<<"("<<min[i].ID1<<","<<min[i].ID2<<","<<min[i].d<<")"<<endl;
        }
        */

         int n,s;
         vector<atom_pair> t1;

           Vector3D Ca;      //calpha atom (DV1)   
           Vector3D GcSc;    //geometrical center of side chain (DV2)
           Vector3D FnAtom;  //First nearest atom to calpha     (DV3)
           Vector3D SnAtom;  //Second nearest atom to calpha    (DV4)
           Vector3D Atom_3rd,Atom_3rd_Ca;
           Vector3D Atom_4th,Atom_4th_Ca;
           Vector3D Atom_5th,Atom_5th_Ca;
           Vector3D Atom_6th,Atom_6th_Ca;
 
           Vector3D FnCa;     //FnAtom-Ca
           Vector3D SnCa;    //SnAtom-Ca
           Vector3D Pn;      // FnCa X SnCa normal vector of plane FnCa X SnCa
           Vector3D SnCa_3rd;
           Vector3D Atom_3rd_4th;
           Vector3D Atom_4th_5th;
           Vector3D Atom_5th_6th;
           Vector3D Atom_6th_Fn;


           Vector3D Tc;       //centroid of triangle from Ca,FnAtom,SnAtom 
           Vector3D TcCa;     //Ca-Tc
           Vector3D ScCa;     //Ca-GCSc
           Vector3D Ave_all_6; 
           double TcScp,TcCaPnp;      //TcCa*ScCa              

        int p1,p2,m,aa_key;    
        double tx,ty,tz,ax,ay,az;  
          
        for(n=1;n<=rows;n++)
        {    
            for(i=min.size();i>=0;i--)
            {
             if(min[i].ID1==n)
             {t1.push_back(min[i]);}
            }
        

   
            for(p1=0;p1<surface_calpha_atom.size();p1++)
           {
              //take calpha atom, geometrical center of main residue
             if(surface_calpha_atom[p1].atomID==n)
             {
              Ca.SetXComponent(surface_calpha_atom[p1].x);
              Ca.SetYComponent(surface_calpha_atom[p1].y);
              Ca.SetZComponent(surface_calpha_atom[p1].z);
           
                
              for(p2=0;p2<w_c.size();p2++)
              {
                  w_c[p2].angle=-100;
             
               if(w_c[p2].res_num==surface_calpha_atom[p1].res_num&&
                  w_c[p2].chain_type==surface_calpha_atom[p1].chain_type&&
                  strcmp(w_c[p2].res_id,surface_calpha_atom[p1].res_id)==0)
                {
                     tx=0;ty=0;tz=0;ax=0;ay=0;az=0;
                     m=0;
                     aa_key=p2; 
                     cout<<w_c[p2].res_id<<endl;
    
                  for(j=0;j<w_c[p2].atom_str.size();j++)
                  {
                     //cout<<w_c[p2].atom_str[j].atom_id<<endl;
                    if(strncmp(w_c[p2].atom_str[j].atom_id,"CA ",3)!=0&&strncmp(w_c[p2].atom_str[j].atom_id,"C  ",3)!=0
                     &&strncmp(w_c[p2].atom_str[j].atom_id,"N  ",3)!=0&&strncmp(w_c[p2].atom_str[j].atom_id,"O  ",3)!=0)
                     {
                       tx=tx+w_c[p2].atom_str[j].x;
                       ty=ty+w_c[p2].atom_str[j].y;
                       tz=tz+w_c[p2].atom_str[j].z;
                       m=m+1;
                      } //end if(strncmp
                  }// end for(j=0

                  if(m!=0)
                  {ax=tx/m;
                   ay=ty/m;
                   az=tz/m;
                   }
                } //end if(w_c[p2]
             }//end for(p2=0..
 
              GcSc.SetXComponent(ax);
              GcSc.SetYComponent(ay);
              GcSc.SetZComponent(az);
            }//end if(surface_calpha..


             //take 1st nearest calpha atom to main alpha atom  
           if(surface_calpha_atom[p1].atomID==t1[0].ID2)
           {
            FnAtom.SetXComponent(surface_calpha_atom[p1].x);
            FnAtom.SetYComponent(surface_calpha_atom[p1].y);
            FnAtom.SetZComponent(surface_calpha_atom[p1].z);              
           }

            //take 2nd nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[1].ID2)
           {
            SnAtom.SetXComponent(surface_calpha_atom[p1].x);
            SnAtom.SetYComponent(surface_calpha_atom[p1].y);
            SnAtom.SetZComponent(surface_calpha_atom[p1].z);              
           }
            
           //take 3rd nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[2].ID2)
           {
            Atom_3rd.SetXComponent(surface_calpha_atom[p1].x);
            Atom_3rd.SetYComponent(surface_calpha_atom[p1].y);
            Atom_3rd.SetZComponent(surface_calpha_atom[p1].z);              
           }

            //take 4th nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[3].ID2)
           {
            Atom_4th.SetXComponent(surface_calpha_atom[p1].x);
            Atom_4th.SetYComponent(surface_calpha_atom[p1].y);
            Atom_4th.SetZComponent(surface_calpha_atom[p1].z);              
           }

            //take 5th nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[4].ID2)
           {
            Atom_5th.SetXComponent(surface_calpha_atom[p1].x);
            Atom_5th.SetYComponent(surface_calpha_atom[p1].y);
            Atom_5th.SetZComponent(surface_calpha_atom[p1].z);              
           }

            //take 6th nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[5].ID2)
           {
            Atom_6th.SetXComponent(surface_calpha_atom[p1].x);
            Atom_6th.SetYComponent(surface_calpha_atom[p1].y);
            Atom_6th.SetZComponent(surface_calpha_atom[p1].z);              
           }

           } //end for(p1=0

             FnCa=FnAtom-Ca;
             SnCa=SnAtom-Ca;
             Atom_3rd_Ca=Atom_3rd-Ca;
             Atom_4th_Ca=Atom_4th-Ca;
             Atom_5th_Ca=Atom_5th-Ca;
             Atom_6th_Ca=Atom_6th-Ca;


           
             Pn=Pn.CrossProduct(FnCa,SnCa);//1stx2nd

             SnCa_3rd=SnCa_3rd.CrossProduct(SnCa,Atom_3rd_Ca);//2ndx3rd

             Atom_3rd_4th=Atom_3rd_4th.CrossProduct(Atom_3rd_Ca,Atom_4th_Ca);//3rdx4th

             Atom_4th_5th=Atom_4th_5th.CrossProduct(Atom_4th_Ca,Atom_5th_Ca);//4thx5th

             Atom_5th_6th=Atom_5th_6th.CrossProduct(Atom_5th_Ca,Atom_6th_Ca);//5thx6th

             Atom_6th_Fn=Atom_6th_Fn.CrossProduct(Atom_6th_Ca,FnCa);//6thx1st


             

                

              /*
             Tc.SetXComponent((Ca.GetXComponent()+FnCa.x+SnCa.x)/3);
             Tc.SetYComponent((Ca.y+FnCa.y+SnCa.y)/3);
             Tc.SetZComponent((Ca.z+FnCa.z+SnCa.z)/3);
               */

             //Tc.SetXComponent((Ca.GetXComponent()+FnAtom.GetXComponent()+SnAtom.GetXComponent())*(1/3));
             //Tc.SetYComponent((Ca.GetYComponent()+FnAtom.GetYComponent()+SnAtom.GetYComponent())*(1/3));
             //Tc.SetZComponent((Ca.GetZComponent()+FnAtom.GetZComponent()+SnAtom.GetZComponent())*(1/3));

              
             Tc=Ca+FnAtom+SnAtom;
             Tc=Tc*0.33333;

             Ave_all_6=Pn+SnCa_3rd+Atom_3rd_4th+Atom_4th_5th+Atom_5th_6th+Atom_6th_Fn;
             Ave_all_6=Ave_all_6*0.16667;



             TcCa=Ca-Tc;
             ScCa=Ca-GcSc;

             //two vector dot product 
             TcScp=TcCa.DotProduct(ScCa);
             TcCaPnp=TcCa.DotProduct(Pn);
             double ScCa_All=ScCa.DotProduct(Ave_all_6);
   
             
             double TcCaL=TcCa.CalculateLength();
             double ScCaL=ScCa.CalculateLength();
             double PnL=Pn.CalculateLength();
             double Ave_all_6L=Ave_all_6.CalculateLength(); 

             double TcScL=TcCaL*ScCaL;
             double TcCaPnL=TcCaL*PnL;
             double ScCa_all_L=ScCaL*Ave_all_6L; 
 
              
            

             double local_angle=180*acos(TcScp/TcScL)/3.14;
             double PnTcCa_angle=180*acos(TcCaPnp/TcCaPnL)/3.14;
             double ScCa_all_angle=180*acos(ScCa_All/ScCa_all_L)/3.14; 
 
             w_c[aa_key].angle=local_angle;
             w_c[aa_key].app_angle=ScCa_all_angle;

             
             sp.push_back(w_c[aa_key]);
 
             
 
             cout<<"Ca= "<<Ca<<endl;
             cout<<"GcSc= "<<GcSc<<endl;
             cout<<"FnAtom= "<<FnAtom<<endl;
             cout<<"SnAtom=" <<SnAtom<<endl;
             cout<<"Atom_3rd= "<<Atom_3rd<<endl;
             cout<<"Atom_4th= "<<Atom_4th<<endl;
             cout<<"Atom_5th= "<<Atom_5th<<endl;
             cout<<"Atom_6th= "<<Atom_6th<<endl;

             /*     
             cout<<"FnCa= "<<FnCa<<endl;
             cout<<"SnCa= "<<SnCa<<endl;
             

             cout<<"Pn= "<<Pn<<endl;
             cout<<"SnCa_3rd= "<<SnCa_3rd<<endl;
             cout<<"Atom_3rd_4th= "<<Atom_3rd_4th<<endl;
             cout<<"Atom_4th_5th= "<<Atom_4th_5th<<endl;
             cout<<"Atom_5th_6th="<<Atom_5th_6th<<endl;
             cout<<"Atom_6th_Fn="<<Atom_6th_Fn<<endl;
             cout<<"Ave_all_6= "<<Ave_all_6<<endl;
             */



             /*
             cout<<"Tc= "<<Tc<<endl;
             cout<<"TcCa= "<<TcCa<<endl;
             cout<<"ScCa= "<<ScCa<<endl;
             cout<<"TcCa*ScCa= "<<TcScp<<endl;
             cout<<"TcCaL= "<<TcCaL<<endl;
             cout<<"ScCaL= "<<ScCaL<<endl;
             cout<<"aa_key= "<<aa_key<<endl;
              */

             //cout<<"local_angle= "<<local_angle<<endl;
             //cout<<setiosflags(ios::left)<<setw(5)<<"PnTcCa_angle= "
               //  <<setiosflags(ios::fixed)
                // <<setiosflags(ios::left)<<setprecision(2)<<setw(10)<<PnTcCa_angle<<endl;

           
         t1.clear();    
        }

/*          
          int neighbors_of_atom;
         
         for(i=0;i<chain_all_atom.size();i++)
         {
          //just count surface atom 
          if(chain_all_atom[i].asa>5.0)
           { 
             neighbors_of_atom=0;
           for(j=0;j<chain_all_atom.size();j++)
            {
             if(distance_matrix[i][j]>0&&distance_matrix[i][j]<=5.0)
             {
              neighbors_of_atom=neighbors_of_atom+1;
             }
            }
          cout<<setw(5)<<chain_all_atom[i].res_num
              <<setw(5)<<chain_all_atom[i].chain_type
              <<setw(5)<<chain_all_atom[i].res_id
              <<setw(5)<<chain_all_atom[i].atom_id
              <<"has "<< neighbors_of_atom << " neighbor atom within 5A "<<endl;
           } // end of if   
             
         }
               


 
          cout<<setw(5)<<chain_all_atom[0].res_num
              <<setw(5)<<chain_all_atom[0].chain_type
              <<setw(5)<<chain_all_atom[0].res_id
              <<setw(5)<<chain_all_atom[0].atom_id
              <<setw(10)<<chain_all_atom[0].x
              <<setw(10)<<chain_all_atom[0].y
              <<setw(10)<<chain_all_atom[0].z
              <<setw(10)<<chain_all_atom[0].asa<<endl;  
        Print(distance_matrix);
         //cout<<distance_matrix[0][0]<<endl;
*/
             
}// end of function res_alpha_atom_

