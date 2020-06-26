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
float omega_f(float,float,float,float,float,float,float,float,float,float,float,float);
float logterm(float,float,float);

#define vert_n 25350
#define face_n 49152
#define vec_n (face_n+4)/2
#define edge_n vec_n+face_n-2

float v[vert_n][3];
int f[face_n][3];
float normals[face_n][3];
float F_fs[face_n][3][3];
float E_es[edge_n][3][3];
float edges_len[edge_n];
float edges[edge_n][3];
float r_vec[vert_n][3];
float r_len[vert_n];


int main()
{
    float vector1[1][3],vector2[1][3];
    ifstream file("itokawa2.txt");
    string content;
    int cont =0;

    while(file >> content) {

      if ((cont%4)!=0 && cont <vert_n*4)
    {v[cont/4][(cont-1)%4]=std::stof(content)*1000;}

       if ( cont >vert_n*4 && (cont-1)%4 !=0)
    {f[(cont-1-4*vert_n)/4][(cont-2)%4]=std::stoi(content);}
      cont++;
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
cont =0;
    for(int i=0;i<face_n;i++)
   {for(int j=0;j<3;j++)
     {
           con0=f[i][j];
           con1=f[i][(j+1)%3];

           common_edge = find_common_edge(con0,con1,i);

           if(common_edge!=-1)

           {float cross1[3],cross2[3];
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
           edges_len[cont]=edge_len;
           edges[cont][0]=con0;edges[cont][1]=con1;edges[cont][2]=i;edges[cont][3]=ii;

            cross1[0] = edge_vec[0][1]*normals[i][2]-edge_vec[0][2]*normals[i][1];
            cross1[1] = edge_vec[0][2]*normals[i][0]-edge_vec[0][0]*normals[i][2];
            cross1[2] = edge_vec[0][0]*normals[i][1]-edge_vec[0][1]*normals[i][0];

            cross2[0] = -edge_vec[0][1]*normals[ii][2]+edge_vec[0][2]*normals[ii][1];
            cross2[1] = -edge_vec[0][2]*normals[ii][0]+edge_vec[0][0]*normals[ii][2];
            cross2[2] = -edge_vec[0][0]*normals[ii][1]+edge_vec[0][1]*normals[ii][0];

              for(int j = 0;j<3;j++)
              { for(int k=0;k<3;k++)
              {
                E_es[cont][j][k]=normals[i][j]*cross1[k]+normals[ii][j]*cross2[k];

              }
               }
             cont ++;
           }

     }

   }
// End of Part 2 -------------------ok ------------------------------//

float field[3]={500.0,0.0,0.0};

for(int i=0;i<vert_n;i++)
{
    for(int j=0;j<3;j++)
    {r_vec[i][j]=v[i][j]-field[j];}

    r_len[i]=norm(r_vec[i][0],r_vec[i][1],r_vec[i][2]);
}

//End of Part 3 ------------------------ok----------------------//

float E_term,F_term=0;
int a,b,c;
float wf;
float prodF[3],prodE[3];

for (int i=0;i<face_n;i++)
{
a=f[i][0]-1;b=f[i][1]-1;c=f[i][2]-1;

wf=omega_f(r_vec[a][0],r_vec[a][1],r_vec[a][2],r_vec[b][0],r_vec[b][1],r_vec[b][2],r_vec[c][0],r_vec[c][1],r_vec[c][2],r_len[a],r_len[b],r_len[c]);

for(int j=0;j<3;j++)
{prodF[j]=0;
for(int k=0;k<3;k++)
  {
    prodF[j] = prodF[j]+F_fs[i][j][k]*r_vec[a][k];
  }
}
F_term=F_term+wf*(r_vec[a][0]*prodF[0]+r_vec[a][1]*prodF[1]+r_vec[a][2]*prodF[2]);

}
//End of Part 4 ------------------------ok----------------------//
float Le;
for (int i=0;i <edge_n;i++)
{edge_len=edges_len[i];
 a=edges[i][0]-1;b=edges[i][1]-1;
 Le=logterm(r_len[a],r_len[b],edge_len);

for(int j=0;j<3;j++)
{prodE[j]=0;
for(int k=0;k<3;k++)
  {
    prodE[j] = prodE[j]+E_es[i][j][k]*r_vec[a][k];
  }
}
E_term=E_term+Le*(r_vec[a][0]*prodE[0]+r_vec[a][1]*prodE[1]+r_vec[a][2]*prodE[2]);
}
//End of Part 5 ------------------------ok----------------------//

float Pot=0.5*1900*(E_term-F_term);
cout<<Pot;
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
float omega_f(float a0,float a1,float a2,float b0,float b1,float b2,float c0,float c1,float c2,float la,float lb,float lc)
{float x,y,wf;

y=a0*(b1*c2-c1*b2)+a1*(b2*c0-c2*b0)+a2*(b0*c1-b1*c0);
x=0;
x=x+la*lb*lc;
x=x+la*(b0*c0+b1*c1+b2*c2)+lb*(a0*c0+a1*c1+a2*c2)+lc*(b0*a0+b1*a1+b2*a2);

wf = 2*atan2(y,x);

}
float logterm(float a,float b,float c)
{float Le= log((a+b+c)/(a+b-c));
return Le;

}
