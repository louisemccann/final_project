(* ::Package:: *)

(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: Louise *)
(* :Date: 2018-04-22 *)

EC0 = 0.214; EJ0 = 52; NG0 = 1*10^-10;

(* Cottet method *)
Em[k0_] := Module[{k=k0}, EC = EC0; EJ=EJ0;NG=NG0;
  				q = -2 EJ/EC;
  				r = k + 1 - Mod[k + 1, 2] + 2 NG (-1)^k;
  				(EC*4) MathieuCharacteristicA[r, q]];

data := Table[Em[x], {x,0,80}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data_ng_nearly0.dat",data ];

NG0=0.0;
data := Table[Em[x], {x,0,80}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data_ng_0.dat",data ];

EC0 = 0.214; EJ0 = 52;
Emkoch[m0_, ng0_] := Module[{m=m0, ng=ng0},EC0; EJ = EJ0;
		k = (Round[Mod[2*ng + (-1/2.0),2]] * (Round[ng] + -1*(-1)^m) *Quotient[m+1,2]) + (Round[Mod[2*ng + (1/2.0),2]] * (Round[ng] + 1*(-1)^m) *Quotient[m+1,2]);
    EC0*MathieuCharacteristicA[(2*(ng+k)),(-EJ0/(2*EC0))]];
wavekoch[m0_,phi0_] := Module[{m=m0, phi=phi0},EC = EC0; EJ = EJ0; NG = NG0;
				q = -EJ/2*EC;
				(Exp[I*NG*phi] / Sqrt[2]) * (MathieuC[4*Em[m]/EC,q,phi/2]+ I (-1)^(m+1) MathieuS[4*Em[m]/EC,q,phi/2])];
				
Plot[Abs[wavekoch[1,x]],{x,-2*Pi,2*Pi}]	
nn=1
NIntegrate[Exp[I*nn*Phi]*wavekoch[1, Phi], {Phi,0,2*Pi}]
Abs[%]^2
N[%]




