function Pot = PolyGrav(f,v,vert_n,face_n,edge_n,F_Fs,E_es,edges_len,edges,x,y,z)

field = [x y z];

r_vec=zeros(vert_n,3);
r_len=zeros(vert_n,1);

for i = 1 : vert_n
    r_vec(i,:)=v(i,:)-field;
    r_len(i)=norm(r_vec(i,:));
end

F_term=0;

for n =1 :face_n
    
    polygon = f(n,:);
    a=polygon(1);
    b= polygon(2);
    c=polygon(3);
    F_f = F_Fs(:,:,i);
%     norma = normal(i,:); 
    
    wf = omega_f([r_vec(a,:);r_vec(b,:);r_vec(c,:)],[r_len(a);r_len(b);r_len(c)]);
    
    F_term=F_term+dot(r_vec(a,:),F_f*r_vec(a,:)')*wf;

end
        
E_term = 0;

for n =1 :edge_n
edge_len = edges_len(n);
    E_e=E_es(:,:,n);
    edge=edges(n,:);
    a=edge(1);b=edge(2);
    E_term =E_term+ dot(r_vec(a,:),E_e*r_vec(a,:)')*logterm(r_len(a),r_len(b),edge_len);
end

Pot = 0.5*1900*(E_term-F_term);

end
 

