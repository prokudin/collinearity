Off Statistics;
.sort;
Index mu,nu;
Vector P,q,Pplus,Pt;
Symbols A,B;
Symbols M,Q,xN;
CFunction F1,F2;
Symbols n;
Local [P.P]  = P(1)*P(2)+P(2)*P(1);
Local [P.q]  = P(1)*q(2)+P(2)*q(1);
id P(1)=Pplus;
id P(2)=M^2/2/Pplus;
id q(1)=-xN*Pplus;
id q(2)=Q^2/2/xN/Pplus;
id Pplus=(Q/sqrt_(2)/xN);
id 1/Pplus=1/(Q/sqrt_(2)/xN);
.sort;
Local [W1]  = -d_(mu,nu)+q(mu)*q(nu)/q.q;
Local [W2]  = (P(mu)-q(mu)*P.q/q.q)*(P(nu)-q(nu)*P.q/q.q)/P.q;
Local [Pg]  = d_(mu,nu);
Local [Ppp] = P(mu)*P(nu);
Local [P]   = A*[Pg] + B*[Ppp];
Local [PW1] = [P]*[W1];
Local [PW2] = [P]*[W2];
id P.P=[P.P];
id q.q^n?=(-Q^2)^n;
id P.q^n?=([P.q])^n;
.sort;
id P.P=[P.P];
id P.q^n?=([P.q])^n;
id q.q^n?=(-Q^2)^n;
.sort;
format fortran;
Print [PW1];
Print [PW2];
.end;
