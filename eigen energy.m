(* Mathematica Source File  *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)
(* :Author: Louise *)
(* :Date: 2018-04-22 *)

EC0 = 0.214; EJ0 = 52; NG0 = 0.000001;

(* Cottet method *)
Em[k0_] := Module[{k=k0}, EC = EC0; EJ=EJ0;NG=NG0;
                  q=-2EJ/EC;
                  r = k + 1 - Mod[k + 1, 2] + 2 NG (-1)^k;
                  (EC*4) MathieuCharacteristicA[r, q]];
data := Table[Em[x]/(2/EJ0), {x,0,80}];
Export["C:\\Users\\Louise\\Documents\\project\\final_project\\data.dat",data ];