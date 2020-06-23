function wf  = omega_f(r,r_len)

y=dot(r(1,:),cross(r(2,:),r(3,:)));
x=0;
x=x+r_len(1)*r_len(2)*r_len(3);
for i = 1 : 3
x=x+r_len(i)*dot(r(mod(i,3)+1,:),r(mod(i+1,3)+1,:));    
end

wf= 2*atan2(y,x);
end