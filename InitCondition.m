tic 
clc,clear all

fname='itokawa2.txt';
d = importdata(fname);
F = 25350;
v =[d(1:F,2) d(1:F,3) d(1:F,4)];
f=[d(F+2:end,2) d(F+2:end,3) d(F+2:end,4)];

 v=v*1000;

vert_n = length(v);
face_n = length(f);
vec_n=(face_n+4)/2;
edge_n = vec_n + face_n - 2;


normal=zeros(face_n,3);
F_Fs=zeros(3,3,face_n);

for i = 1:face_n
vector1=v(f(i,2),:)-v(f(i,1),:);
vector2=v(f(i,3),:)-v(f(i,2),:);
% edge_len=norm(vector1(i,:));
 vector1=vector1/norm(vector1);
 vector2=vector2/norm(vector2);

normal(i,:) = cross(vector1,vector2);
normal(i,:)=normal(i,:)/norm(normal(i,:));

 F_Fs(:,:,i) = normal(i,:)'*normal(i,:);

end

% common_edge=zeros(face_n,1);

 E_es=zeros(3,3,edge_n);
 edges_len=zeros(edge_n,1);
 edges=zeros(edge_n,4);

count=0;
for i = 1: face_n
  for j = 1:3
 con(1) = f(i,j);
 con(2) = f(i,mod(j,3)+1);

 common_edge= find_common_edge(f,face_n,con,i);

 if common_edge ~= -1
     
     count=count+1;
     ii = common_edge;
     p1 = v(con(1),:);
     p2 = v(con(2),:);
     edge_vec = p2 - p1;
     edge_len = norm(edge_vec);
     edge_vec = edge_vec / edge_len;
 
     E_es(:,:,count)= normal(i,:)'*cross(edge_vec,normal(i,:))+normal(ii,:)'*cross(-edge_vec,normal(ii,:));
     edges_len(count)=edge_len;
     edges(count,:)=[con(1) con(2) i ii];
     
 end
  end
end

N =1000;

rmin = 500;
rmax=700;h= 400;

count =1;

while count<=N
 
x(count) = 2*rand*h-h;
z(count) = 2*rmax*rand-rmax;
y(count) = 2*rmax*rand-rmax;
r = sqrt(y(count)^2+z(count)^2);

count = count+1;
%  
if r<rmin|r>rmax
    count = count - 1;
     continue;
end

end

R=zeros(1,N);
Pot=zeros(1,N);
K=zeros(1,N);
G= 1;
M = 3.5e10;

for i = 1:N
    
R(i)= sqrt(x(i)^2+y(i)^2+z(i)^2); 
Pot(i) = PolyGrav(f,v,vert_n,face_n,edge_n,F_Fs,E_es,edges_len,edges,x(i),y(i),z(i));
K(i)=(G*M/R(i));

end


plot(R,Pot,'*')
hold on
plot(R,K,'k.')





%  plot3(x,y,z,'.')
toc