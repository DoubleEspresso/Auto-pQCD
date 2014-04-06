off statistics;
autodeclare index dum,tmp;
Indices mu,p,k;
CFunctions p,k,kh,sin,cos,M,S,B,b,d,SUM,sum,I,F,f,sig,eps,s,m;
CFunctions J1,...,J1000;
Function G,g,T,t;
Symbol g0,AA,aa,rh,BB,bb,rw,m0,DB,db,DQ,dq,a,pow,sf;
Local str = (-i_)*cos((a*(Subscript(k,mu)+Subscript(p,mu)))/2)*G(mu)-rw*sin((a*(Subscript(k,mu)+Subscript(p,mu)))/2);
repeat;
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?)*G(dum8?) = G(dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8);
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?)*G(dum7?) = G(dum1,dum2,dum3,dum4,dum5,dum6,dum7);
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?)*G(dum6?) = G(dum1,dum2,dum3,dum4,dum5,dum6);
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?)*G(dum5?) = G(dum1,dum2,dum3,dum4,dum5);
id G(dum1?)*G(dum2?)*G(dum3?)*G(dum4?) = G(dum1,dum2,dum3,dum4);
id G(dum1?)*G(dum2?)*G(dum3?) = G(dum1,dum2,dum3);
id G(dum1?)*G(dum2?) = G(dum1,dum2);
endrepeat;
contract;
print str;
*********
*ANSWER
.end
