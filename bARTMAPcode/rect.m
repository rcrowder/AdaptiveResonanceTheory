%function A=rect(B,t), rectifies to value above t.
function A=rect(B,t);
A=(B-t+abs(B-t))/2;
A=A+sign(A).*t;
