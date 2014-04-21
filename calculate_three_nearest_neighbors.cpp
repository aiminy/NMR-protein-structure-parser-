/*
this program calculate three nearest calpha neighboring atoms to seed calpha atom
and construct a triangle by using three nearest calpha neighboring atoms
*/


#include "cal_tot_asa_of_res.h"
#include "tmatrix.h"
#include "Vector3D.h"
//#include "inter_atom_distance_matrix.cpp"



template <class T>


int partitioni1(vector<double> &, int, int);
void quicksort1(vector<double> &, int, int);

int partition1(vector <atom_pair> &array, int top, int bottom)
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

void quicksort1(vector <atom_pair> &array, int top, int bottom)
{
     int middle;
     if (top < bottom)
    {
          middle = partition1(array, top, bottom);
          quicksort1(array, top, middle);   // sort top partition
          quicksort1(array, middle+1, bottom);    // sort bottom partition
     }
     return;
}


void calculate_three_nearest_neighbors(vector<aa> &w_c,vector<aa> &sp)
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

 
        quicksort1(min, 0, min.size());
        
         
        //cout<<"After sorting"<<endl;
        //for(i=0;i<min.size();i++)
        //{
        // cout<<"("<<min[i].ID1<<","<<min[i].ID2<<","<<min[i].d<<")"<<endl;
        //}
       

         int n,s;
         vector<atom_pair> t1;

           Vector3D Ca;       //the calpha atom of seed residue   
           Vector3D GcSc;     //the geometrical center of side chain of seed residue 
           Vector3D FnAtom;   //the first nearest calpha atom to calpha atom of seed residue
           Vector3D SnAtom;   //the second nearest calpha atom to calpha atom of seed residue
           Vector3D Atom_3rd; //the third nearest calpha atom to calpha atom of seed residue


           Vector3D Sn_Fn;       //SnAtom-FnAtom
           Vector3D Atom_3rd_Fn; //Atom_3rd-FnAtom

           Vector3D Tc,TcCa,ScCa; 
           Vector3D Pn;      // normal vector of triangle that contains FnAtom,SnAtom and Atom_3rd

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
                  w_c[p2].angle_atom_3c=-100;
             
               if(w_c[p2].res_num==surface_calpha_atom[p1].res_num&&
                  w_c[p2].chain_type==surface_calpha_atom[p1].chain_type&&
                  strcmp(w_c[p2].res_id,surface_calpha_atom[p1].res_id)==0)
                {
                     tx=0;ty=0;tz=0;ax=0;ay=0;az=0;
                     m=0;
                     aa_key=p2; 
                     //cout<<w_c[p2].res_id<<endl;
    
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



             Sn_Fn=SnAtom-FnAtom;
             Atom_3rd_Fn=Atom_3rd-FnAtom;

             Pn=Pn.CrossProduct(Sn_Fn,Atom_3rd_Fn);//Sn_Fn x Atom_3rd_Fn

              
             Tc=FnAtom+SnAtom+Atom_3rd;
             Tc=Tc*0.33333;


             TcCa=Ca-Tc;
             ScCa=Ca-GcSc;

             //two vector dot product 
             TcScp=Pn.DotProduct(ScCa);

             double ScCaL=ScCa.CalculateLength();
             double PnL=Pn.CalculateLength();

             double TcScL=PnL*ScCaL;
             double PnTcCa_angle=180*acos(TcScp/TcScL)/3.14;

             w_c[aa_key].angle_atom_3c=PnTcCa_angle;
             
             sp.push_back(w_c[aa_key]);
 
         t1.clear();    
        }
      }

     }// end of function calculate_three_nearest_neighbors

