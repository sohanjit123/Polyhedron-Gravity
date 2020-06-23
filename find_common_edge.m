function out = find_common_edge(f,face_n,con,i)

 n = i + 1;
    
 for i = n :face_n
     for j = 1 : 3
        
         if con(2) == f(i,j) && con(1) == f(i,mod(j,3)+1)
             out=i; return
         end
 
     end  
 end
 out = -1;
end
