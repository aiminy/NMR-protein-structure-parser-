#include "cal_tot_asa_of_res.h"
#include "tmatrix.h"
#include "Vector3D.h"
//#include "triangle.h"
#include "linked_list.h"

//template <class T>

int partition(vector<double> &, int, int);
void quicksort(vector<double> &, int, int);
double angle_between_two_vectors(Vector3D A,Vector3D B);

double angle_between_two_vectors(Vector3D A,Vector3D B)
{
  double A_DOT_B=A.DotProduct(B);
  double A_length=A.CalculateLength();
  double B_length=B.CalculateLength();
  double AB=A_length*B_length;
  double angle=180*acos(A_DOT_B/AB)/3.14;
  return angle;
}

int partition(vector <atom_pair> &array, int top, int bottom)
{
     atom_pair x1;

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
        int i,j,rows,cols,num_total_atom,num_of_calpha,rows1,cols1;
            

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


        //center of structure

          double ws_x,ws_y,ws_z;
          double how_many_atom;      
          ws_x=0; 
          ws_y=0; 
          ws_z=0; 
          how_many_atom=0;

       for(i=0;i<w_c.size();i++)
         {
                  w_c[i].global_angle=0;  
            
          for(j=0;j<w_c[i].atom_str.size();j++)
           {
            ws_x=ws_x+w_c[i].atom_str[j].x;
            ws_y=ws_y+w_c[i].atom_str[j].y;
            ws_z=ws_z+w_c[i].atom_str[j].z;
            how_many_atom=how_many_atom+1;
           }
         }

            ws_x=ws_x/how_many_atom;
            ws_y=ws_y/how_many_atom;
            ws_z=ws_z/how_many_atom;

         //cout<<ws_x<<" "<<ws_y<<" "<<ws_z<<" "<<endl;   

         vector<atom> surface_calpha_atom;
         int k=1;
         for(i=0;i<chain_all_atom.size();i++)
          {
           //if((strcmp(chain_all_atom[i].atom_id,"CA ")==0)&&chain_all_atom[i].asa>5&&
           if(strncmp(chain_all_atom[i].res_id,"GLY",3)!=0)
             {
              chain_all_atom[i].atomID=k; 
              surface_calpha_atom.push_back(chain_all_atom[i]);
              k++;
             }//end if
          }//end for

        
           
	//consider all surface atom
         vector<atom> surface_atom;  
         int k1=1;
 
         for(i=0;i<chain_all_atom.size();i++)
          {
           //if(chain_all_atom[i].asa>5
            if(strcmp(chain_all_atom[i].atom_id,"CA ")==0)
             //&&strncmp(chain_all_atom[i].res_id,"GLY",3)!=0)
             {
              chain_all_atom[i].atomID=k1; 
              surface_atom.push_back(chain_all_atom[i]);
              k1++;
             }//end if
          }//end for

         rows1=cols1=surface_atom.size();
        
         tmatrix<atom_pair> distance_matrix1(rows1,cols1);
     
         double t_all_surface_atom;
         
         t_all_surface_atom=0; 
         for(i=0;i<surface_atom.size();i++)
         {
           for(j=0;j<surface_atom.size();j++)
            {
             distance_matrix1[i][j].ID1=surface_atom[i].atomID;
             distance_matrix1[i][j].ID2=surface_atom[j].atomID;

             distance_matrix1[i][j].d=sqrt(pow((surface_atom[i].x-surface_atom[j].x),2)+
                    pow((surface_atom[i].y-surface_atom[j].y),2)+
                    pow((surface_atom[i].z-surface_atom[j].z),2));
            }
         }

        
         vector<atom_pair> min1;
         
        for(i=0;i<rows1;i++)
        {
         for(j=0;j<cols1;j++)
          {
           if(distance_matrix1[i][j].d!=0)
           {
           min1.push_back(distance_matrix1[i][j]);
          }
          }
       } 

        quicksort(min1, 0, min1.size());

         int n_all_surface_atom,s_all_surface_atom;
         vector<atom_pair> t1_all_surface_atom;
        
        for(n_all_surface_atom=1;n_all_surface_atom<=rows1;n_all_surface_atom++)
        {    
            for(i=min1.size();i>=0;i--)
            {
             if(min1[i].ID1==n_all_surface_atom)
             {t1_all_surface_atom.push_back(min1[i]);}
            }
       }
        
        //cout<<"After sorting"<<endl;
        //for(i=0;i<t1_all_surface_atom.size();i++)
        //{
        //if(t1_all_surface_atom[i].d<7.0)
  
       //{cout<<"("<<t1_all_surface_atom[i].ID1<<","<<t1_all_surface_atom[i].ID2<<","<<t1_all_surface_atom[i].d<<")"<<endl;
       // }

       // }
        


          
         rows=cols=surface_calpha_atom.size();
        
         tmatrix<atom_pair> distance_matrix(rows,cols);
     
         double t;
         
         t=0; 
         for(i=0;i<surface_calpha_atom.size();i++)
         {
           for(j=0;j<surface_calpha_atom.size();j++)
            {
             distance_matrix[i][j].ID1=surface_calpha_atom[i].atomID;
             distance_matrix[i][j].ID2=surface_calpha_atom[j].atomID;

             //cout<<distance_matrix[i][j].ID1<<" "<<distance_matrix[i][j].ID2<<endl; 

             distance_matrix[i][j].d=sqrt(pow((surface_calpha_atom[i].x-surface_calpha_atom[j].x),2)+
                    pow((surface_calpha_atom[i].y-surface_calpha_atom[j].y),2)+
                    pow((surface_calpha_atom[i].z-surface_calpha_atom[j].z),2));
            }
         }

        
         vector<atom_pair> min;
         
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

 
        quicksort(min, 0, min.size());
        
            //for(i=min.size();i>=0;i--)
            //{
            // if(min[i].ID1==n)
             //{t1.push_back(min[i]);}
              //cout<<" ID1 "<<min[i].ID1<<" ID2 "<<min[i].ID2<<" d "<<min[i].d<<endl;  
            //}
         

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

           // triangle by FnAtom,SnAtom,Atom_3rd

           Vector3D Sn_Fn;       // SnAtom-FnAtom
           Vector3D Atom_3rd_Fn; // Atom_3rd-FnAtom
           Vector3D Pn_3;        // normal vector of triangle by FnAtom,SnAtom,Atom_3rd
           Vector3D Pn_3_re;
           Vector3D Tc_3;        // centroid of triangle by FnAtom,SnAtom,Atom_3rd  
           Vector3D Tc_3_Ca;     // Ca-Tc_3


           Vector3D Tc;       //centroid of triangle from Ca,FnAtom,SnAtom 
           Vector3D TcCa;     //Ca-Tc
           Vector3D ScCa;     //Ca-GCSc
           Vector3D Ave_all_6; 
           Vector3D Plane_Normal;
           Vector3D triangle_center;
           Vector3D t_add_n;
           Vector3D w_t_c;

           Point v1,v2,v3;

           Triangle triangle;
           node<Triangle> *head_node,*previous_node;
           previous_node=NULL;
           
                      
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
        
         //cout<<n<<endl;    
	for(p1=0;p1<surface_calpha_atom.size();p1++)
           {
              //take calpha atom, geometrical center of main residue
               //Vector3D Ca;
             if(surface_calpha_atom[p1].atomID==n)
             {
              //cout<<surface_calpha_atom[p1].x<<endl;   
              Ca.SetXComponent(surface_calpha_atom[p1].x);
              Ca.SetYComponent(surface_calpha_atom[p1].y);
              Ca.SetZComponent(surface_calpha_atom[p1].z);
           
              //cout<<"Ca"<<Ca<<endl;
   
              for(p2=0;p2<w_c.size();p2++)
              {
                  w_c[p2].angle=-100;
                  //w_c[p2].global_angle=0;  
             
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
              //cout<<"GcSc"<<GcSc<<endl;
 
            }//end if(surface_calpha..



             //take 1st nearest calpha atom to main alpha atom  
           if(surface_calpha_atom[p1].atomID==t1[0].ID2)
           {
            FnAtom.SetXComponent(surface_calpha_atom[p1].x);
            FnAtom.SetYComponent(surface_calpha_atom[p1].y);
            FnAtom.SetZComponent(surface_calpha_atom[p1].z);
            //cout<<t1[0].ID2<<endl;
            v1.SetPointIndex(t1[0].ID2);           
           }

            //take 2nd nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[1].ID2)
           {
            SnAtom.SetXComponent(surface_calpha_atom[p1].x);
            SnAtom.SetYComponent(surface_calpha_atom[p1].y);
            SnAtom.SetZComponent(surface_calpha_atom[p1].z);
            //cout<<t1[1].ID2<<endl;
            v2.SetPointIndex(t1[1].ID2);           
           }
            
           //take 3rd nearest calpha atom to main alpha atom           
           if(surface_calpha_atom[p1].atomID==t1[2].ID2)
           {
            Atom_3rd.SetXComponent(surface_calpha_atom[p1].x);
            Atom_3rd.SetYComponent(surface_calpha_atom[p1].y);
            Atom_3rd.SetZComponent(surface_calpha_atom[p1].z);
            //cout<<t1[2].ID2<<endl;
            v3.SetPointIndex(t1[2].ID2);           
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


             //vector point from gcsc to ca  
             ScCa=Ca-GcSc;

             double x1,y1,z1,x2,y2,z2,x3,y3,z3,A,B,C,D;
             x1=FnAtom.GetXComponent();
             y1=FnAtom.GetYComponent();
             z1=FnAtom.GetZComponent();
                 
             v1.SetPointCoordinate(x1,y1,z1);
             //triangle.first.x=x1;
             //triangle.first.y=y1;
             //triangle.first.z=z1;
             
 
             x2=SnAtom.GetXComponent();
             y2=SnAtom.GetYComponent();
             z2=SnAtom.GetZComponent();

             v2.SetPointCoordinate(x2,y2,z2);

             //triangle.Second.x=x2;
             //triangle.second.y=y2;
             //triangle.Second.z=z2;
              

             x3=Atom_3rd.GetXComponent();
             y3=Atom_3rd.GetYComponent();
             z3=Atom_3rd.GetZComponent();

             v3.SetPointCoordinate(x3,y3,z2);
             //v3.SetPointIndex(12);           

             //cout<<v1;
             //cout<<v2;
             //cout<<v3; 
 
             //cout<<"here is triangle"<<endl;

             triangle.SetTriangleThreePoint(v1,v2,v3);
                
             //cout<<triangle<<endl;

             //insert_an_node(head_node,previous_node,triangle);

             //cout<<"head"<<head_node->data<<endl;
             
             //cout<<"previous"<<previous_node->data<<endl;

             //triangle.third.x=x3;
             //triangle.third.y=y3;
             //triangle.third.z=z3;


             //triangle center by FnAtom,SnAtom,Atom_3rd
             double t_c_x,t_c_y,t_c_z;
 
             t_c_x=(x1+x2+x3)/3;
             t_c_y=(y1+y2+y3)/3;
             t_c_z=(z1+z2+z3)/3;

             //create 3D Vector: whole structure center -> t_c_of_triangle
             double w_t_c_x,w_t_c_y,w_t_c_z;

             w_t_c_x=t_c_x-ws_x;
             w_t_c_y=t_c_y-ws_y;
             w_t_c_z=t_c_z-ws_z;

              
             w_t_c.SetXComponent(w_t_c_x);
             w_t_c.SetYComponent(w_t_c_y);
             w_t_c.SetZComponent(w_t_c_z);
     
            
             
             //create a artificial atom    
	     atom t_c;
             t_c.atom_id =new char[5];
             t_c.atom_id="T_C";
             t_c.x=t_c_x;
             t_c.y=t_c_y;
             t_c.z=t_c_z;
             t_c.asa=0;
             t_c.b=0; 
             w_c[aa_key].atom_str.push_back(t_c);
             

             triangle_center.SetXComponent(t_c_x);
             triangle_center.SetYComponent(t_c_x);
             triangle_center.SetZComponent(t_c_x);

             //triangle_center.Normalize();
  

             A=y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2);
             B=z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2);
             C=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
             D=(-1)*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1)); 
  
             Plane_Normal.SetXComponent(A);
             Plane_Normal.SetYComponent(B);
             Plane_Normal.SetZComponent(C);
             
	     //angle between w_t_c and Plane_Normal
             double global_angle=angle_between_two_vectors(w_t_c,Plane_Normal);
              
             //double length_plane_normal;
             //length_plane_normal=Plane_Normal.CalculateLength();
 
             //Plane_Normal.Normalize();     

             t_add_n=triangle_center+Plane_Normal;

              
            // t_add_n.Normalize();
    
             double t_add_n_x,t_add_n_y,t_add_n_z;
              
             t_add_n_x=t_add_n.GetXComponent();
             t_add_n_y=t_add_n.GetYComponent();
             t_add_n_z=t_add_n.GetZComponent();
                                        
             //create another artificial atom    
	     atom n_t;
             n_t.atom_id =new char[5];
             n_t.atom_id="N_T";
             n_t.x=t_add_n_x;
             n_t.y=t_add_n_y;
             n_t.z=t_add_n_z;
             n_t.asa=0;
             n_t.b=0; 
             w_c[aa_key].atom_str.push_back(n_t);
              
             double x,y,z,sign_calpha,sign_gcsc,x_gcsc,y_gcsc,z_gcsc;
             x=Ca.GetXComponent();
             y=Ca.GetYComponent();
             z=Ca.GetZComponent();

             x_gcsc=GcSc.GetXComponent();
             y_gcsc=GcSc.GetYComponent();
             z_gcsc=GcSc.GetZComponent();
             

             sign_calpha=A*x+B*y+C*z+D;
             sign_gcsc=A*x_gcsc+B*y_gcsc+C*z_gcsc+D;

     
             //cout<<"sign_calpha"<<sign_calpha<<"sign_gcsc"<<sign_gcsc<<endl;  

            
             Pn=Pn.CrossProduct(FnCa,SnCa);//1stx2nd

             SnCa_3rd=SnCa_3rd.CrossProduct(SnCa,Atom_3rd_Ca);//2ndx3rd

             Atom_3rd_4th=Atom_3rd_4th.CrossProduct(Atom_3rd_Ca,Atom_4th_Ca);//3rdx4th

             Atom_4th_5th=Atom_4th_5th.CrossProduct(Atom_4th_Ca,Atom_5th_Ca);//4thx5th

             Atom_5th_6th=Atom_5th_6th.CrossProduct(Atom_5th_Ca,Atom_6th_Ca);//5thx6th

             Atom_6th_Fn=Atom_6th_Fn.CrossProduct(Atom_6th_Ca,FnCa);//6thx1st

              
             Tc=Ca+FnAtom+SnAtom;
             Tc=Tc*0.33333;

             Ave_all_6=Pn+SnCa_3rd+Atom_3rd_4th+Atom_4th_5th+Atom_5th_6th+Atom_6th_Fn;
             Ave_all_6=Ave_all_6*0.16667;



             TcCa=Ca-Tc;

            
             //angle between TcCa and ScCa
             double local_angle=angle_between_two_vectors(TcCa,ScCa);

             //angle between TcCa and Pn 
             double PnTcCa_angle=angle_between_two_vectors(TcCa,Pn);

             //angle between ScCa and Ave_all_6
             double ScCa_all_angle=angle_between_two_vectors(ScCa,Ave_all_6);
             
             FnCa=FnAtom-Ca;
             SnCa=SnAtom-Ca;
             Atom_3rd_Ca=Atom_3rd-Ca;
             Atom_4th_Ca=Atom_4th-Ca;
             Atom_5th_Ca=Atom_5th-Ca;
             Atom_6th_Ca=Atom_6th-Ca;


             //triangle from FnAtom,SnAtom,Atom_3rd;

             Sn_Fn=SnAtom-FnAtom;
             Atom_3rd_Fn=Atom_3rd-FnAtom;

             Pn_3=Pn_3.CrossProduct(Sn_Fn,Atom_3rd_Fn);
             Pn_3_re=Pn_3_re.CrossProduct(Atom_3rd_Fn,Sn_Fn);

             double distance1=abs((A*x+B*y+C*z+D)/sqrt(A*A+B*B+C*C));
             double distance2=abs((A*x_gcsc+B*y_gcsc+C*z_gcsc+D)/sqrt(A*A+B*B+C*C));

             Tc_3=FnAtom+SnAtom+Atom_3rd;
             Tc_3=Tc_3*0.33333;
             Tc_3_Ca=Ca-Tc_3;

             
             //angle between Tc_3_Ca and ScCa
             double angle_3=angle_between_two_vectors(Tc_3_Ca,ScCa);
             //angle between Tc_3_Ca and Pn_3
             double angle_4=angle_between_two_vectors(Tc_3_Ca,Pn_3);

             //angle between ScCa and Pn_3
             double angle_5;
 
             if(global_angle>=90)
             {Pn_3=Pn_3*(-1);  
             angle_5=angle_between_two_vectors(ScCa,Pn_3);
             }
             else
             {
             angle_5=angle_between_two_vectors(ScCa,Pn_3);
             }
      


             //angle between ScCa and Pn_3_re  
             double angle_6=angle_between_two_vectors(ScCa,Pn_3_re);
             double angle_7=180-angle_5; 
                                
             //cout<<w_c[aa_key].res_id<<" "<<angle_3<<" "<<angle_4<<" "<<angle_5<<" "<<angle_6<<" "<<angle_7<<endl;
           
             double angle_8,angle_9,angle_10,angle_11;

             if(sign_calpha>0&&sign_gcsc>0)
             {
             angle_8=angle_between_two_vectors(ScCa,Plane_Normal);
             }
            
             if(sign_calpha<0&&sign_gcsc<0)
             {
             angle_9=angle_between_two_vectors(ScCa,Plane_Normal);

             }
               
             if(sign_calpha>0&&sign_gcsc<0)
             {
             angle_10=angle_between_two_vectors(ScCa,Plane_Normal);

             }

             if(sign_calpha<0&&sign_gcsc>0)
             {
             angle_11=angle_between_two_vectors(ScCa,Plane_Normal);
             }


             w_c[aa_key].global_angle=global_angle;  
             w_c[aa_key].angle=local_angle;
             w_c[aa_key].app_angle=ScCa_all_angle;
             w_c[aa_key].angle_atom_3c=angle_3;
             w_c[aa_key].angle_atom_3n=angle_5;
             w_c[aa_key].distance_calpha_to_3atom_plane=distance1;
             w_c[aa_key].distance_scgc_to_3atom_plane=distance2;
             //w_c[aa_key].distance_scgc_to_3atom_plane=distance2;
             //w_c[aa_key].angle_8=angle_8;
             //w_c[aa_key].angle_9=angle_9;
             //w_c[aa_key].angle_10=angle_10;
             //w_c[aa_key].angle_11=angle_11;
             w_c[aa_key].sign_calpha=sign_calpha;
             w_c[aa_key].sign_gcsc=sign_gcsc;
               //cout<<"CA"<<w_c[aa_key].Ca<<endl;
 

             //cout<<" distance1 is "<<distance1<<"distance2 is "<<distance2<<"angle_5 is"<<angle_5<<endl;
             
             sp.push_back(w_c[aa_key]);
             cout<<"Ca= "<<Ca<<endl;
             cout<<"GcSc= "<<GcSc<<endl;
             cout<<"FnAtom= "<<FnAtom<<endl;
             cout<<"SnAtom=" <<SnAtom<<endl;
             cout<<"Atom_3rd= "<<Atom_3rd<<endl;
             cout<<"Atom_4th= "<<Atom_4th<<endl;
             cout<<"Atom_5th= "<<Atom_5th<<endl;
             cout<<"Atom_6th= "<<Atom_6th<<endl;

           
         t1.clear();    
        }



         //display_node_list(head_node);    

}// end of function res_alpha_atom_
