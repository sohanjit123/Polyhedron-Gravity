#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/array.hpp>
#include "poly_grav.hpp"
#include <math.h>

using namespace std;
float norm(float,float,float);
int find_common_edge(int,int,int);

#define vert_n 25350
#define face_n 49152
#define vec_n (face_n+4)/2
#define edge_n vec_n+face_n-2

float v[vert_n][3];
float f[face_n][3];
float normals[face_n][3];
float F_fs[face_n][3][3];
float E_es[edge_n][3][3];
float edges_len[edge_n];
float edges[edge_n][3];


int main()
{
    float vector1[1][3],vector2[1][3];
    ifstream file("itokawa2.txt");
    string content;
    int count =0;

    while(file >> content) {
    //cout << content << ' ';

    if ((count%4)!=0 && count <vert_n*4)
    {v[count/4][(count-1)%4]=std::stof(content)*1000;}

     if ( count >vert_n*4 && (count-1)%4 !=0)
    {f[(count-1-4*vert_n)/4][(count-2)%4]=std::stod(content);}
      count++;
    }

// End of Part 0------------------OK---------------------//

for(int i = 0;i<face_n;i++)
  {int a=f[i][0]-1 ; int b= f[i][1]-1;int c= f[i][2]-1;
   for(int j=0;j<3;j++)
   {
    vector1[0][j]=v[b][j]-v[a][j];
    vector2[0][j]=v[c][j]-v[b][j];
   }
   float n1 = norm(vector1[0][0],vector1[0][1],vector1[0][2]);
   float n2 = norm(vector2[0][0],vector2[0][1],vector2[0][2]);

   for(int j=0;j<3;j++)
   {
     vector1[0][j]=vector1[0][j]/n1;
     vector2[0][j]=vector2[0][j]/n2;
   }

   normals[i][0]=(vector1[0][1]*vector2[0][2]-vector1[0][2]*vector2[0][1]);
   normals[i][1]=(vector1[0][2]*vector2[0][0]-vector1[0][0]*vector2[0][2]);
   normals[i][2]=(vector1[0][0]*vector2[0][1]-vector1[0][1]*vector2[0][0]);

   float n = norm(normals[i][0],normals[i][1],normals[i][2]);


   for(int j=0;j<3;j++)
   {
     normals[i][j]=normals[i][j]/n;
   }

   for(int j = 0;j<3;j++)
   { for(int k=0;k<3;k++)
      {
        F_fs[i][j][k]=normals[i][j]*normals[i][k];

      }
   }
  }
//End of Part 1-------------------OK--------------------------//
 int con0,con1,common_edge,ii;
 float p1[1][3],p2[1][3],edge_vec[1][3];
 float edge_len;
count =0;
    for(int i=0;i<face_n;i++)
   {for(int j=0;j<3;j++)
     {
           con0=f[i][j];
           con1=f[i][(j+1)%3];

           common_edge = find_common_edge(con0,con1,i);

           if(common_edge!=-1)

           {
               ii = common_edge;

           for(int j=0;j<3;j++)
           { p1[0][j]=v[con0-1][j];
             p2[0][j]=v[con1-1][j];
             edge_vec[0][j]=p2[0][j]-p1[0][j];
           }

           edge_len = norm(edge_vec[0][0],edge_vec[0][1],edge_vec[0][2]);

                for(int j=0;j<3;j++)
           { edge_vec[0][j]=edge_vec[0][j]/edge_len;

           }
           edges_len[count]=edge_len;
           edges[count][0]=con0;edges[count][1]=con1;edges[count][2]=i;edges[count][3]=ii;

             count ++;
           }

     }

   }

    return 0;
}








float norm(float a,float b,float c)
{
  float n= sqrt(pow(a,2)+pow(b,2)+pow(c,2));
  return n;
}
int find_common_edge(int con0,int con1,int i)
{
   int n=i+1;

   for(int i=n;i<face_n;i++)
   { for(int j=0;j<3;j++)
      {
         if (con0==f[i][(j+1)%3] &&  con1==f[i][j])
          {  return i;}

      }
    }
return -1;
}
