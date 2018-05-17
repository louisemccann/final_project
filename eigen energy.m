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
    
data2 := Table[Emkoch[x,0.0],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data2.dat",data2 ];
    
data2 := Table[Emkoch[x,0],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_0.dat",data2 ];
data2 := Table[Emkoch[x,0.1],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_01.dat",data2 ];
data2 := Table[Emkoch[x,0.2],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_02.dat",data2 ];
data2 := Table[Emkoch[x,0.3],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_03.dat",data2 ];
data2 := Table[Emkoch[x,0.4],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_04.dat",data2 ];
data2 := Table[Emkoch[x,0.5],{x,0,40}];
{
 {Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_05.dat",data2 ];},
 {data2 := Table[Emkoch[x,0.6],{x,0,40}];
  Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data\\data_ng_06.dat",data2 ];},
 {\[Placeholder]}
}

wavekoch[m0_] := Module[{m=m0},EC = EC0; EJ = EJ0; NG = ng0;
				r = -2 * (NG- Emkoch[m0])
				q = -EJ/2*EC
				Table[Module[{},(Exp[I*NG*phi] / Sqrt[2]) * MathieuC[r,q,phi/2]], {phi, -2*Pi, 2*Pi}]];

EC0 = 0.214; EJ0 = 52; ng0 = 1E-5;
Emkoch[m0_] := Module[{m=m0},EC0; EJ = EJ0; NG = ng0;
		k = (Round[Mod[2*NG + -1/2.0,2]] * (Round[NG] + -1*(-1)^m) *Quotient[m+1,2]) + (Round[Mod[2*NG + 1/2.0,2]] * (Round[NG] + 1*(-1)^m) *Quotient[m+1,2]);
    EC0*MathieuCharacteristicA[(2*NG+k),(-EJ0/(2*EC0))]];
data2 := Table[Emkoch[x],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data2.dat",data2 ];


