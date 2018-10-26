--we choose the rings
G=frac(QQ[Z]/(Z^2-2641))
R=G[X_0,X_1,X_2]
--we constract the four linear forms L_1,..,L_4
(a,b,c)=((Z+119)/(4*(Z+47)),(-3*(Z+39))/16,1)
(d,e,f)=(-2*(Z+59)/(Z+47),5,1)
(g,h,i)=(2,9,1)
(j,k,l)=(1,(Z+11)/(4),1)
L1=a*X_0+b*X_1+c*X_2
L2=d*X_0+e*X_1+f*X_2
L3=g*X_0+h*X_1+i*X_2
L4=j*X_0+k*X_1+l*X_2
SL={L1,L2,L3,L4}
--we constract the star configuration X(L_1,...,L_4)
IX=ideal for t from 1 to 4 list product toList(set(SL)-set{SL_(t-1)})
M=X_0*X_1*X_2
--we find the perp ideal of M
Mperp=sub(inverseSystem(sub(M,QQ[X_0,X_1,X_2])),R)
use R
--it shwos that I_X(4) is in Mperp
isSubset(IX,Mperp)
--we find the matrix of the coef of the four linear forms L_1,..,L_4
Matcoefs=transpose matrix{{L1},{L2},{L3},{L4}}
(W,C) = coefficients(Matcoefs, Monomials=>{X_0,X_1,X_2})
--here we check that ther are linearly independent
minors(3,transpose sub(C,G))
--here we verify that the weak Hadamard star configuration X(L_1,...,L_4) is a Hadamard star configuration.!
rank transpose matrix for a in entries sub(C,G) list apply(a, i-> 1/i)
