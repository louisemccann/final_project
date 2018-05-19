(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: Louise *)
(* :Date: 2018-05-19 *)

EC0 = 0.214; EJ0 = 52; ng0 = 0;
Emkoch[m0_] := Module[{m=m0},EC0; EJ = EJ0; NG = ng0;
		k = (Round[Mod[2*NG + (-1/2.0),2]] * (Round[NG] + -1*(-1)^m) *Quotient[m+1,2]) + (Round[Mod[2*NG + (1/2.0),2]] * (Round[NG] + 1*(-1)^m) *Quotient[m+1,2]);
    EC0*MathieuCharacteristicA[(2*NG+k),(-EJ0/(2*EC0))]];
data2 := Table[Emkoch[x],{x,0,40}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data2.dat",data2 ];