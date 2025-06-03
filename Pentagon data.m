(* ::Package:: *)

(* ::Chapter:: *)
(*5-point letters*)


Quit[]


(* ::Section::Initialization::Closed:: *)
(*Two-loop planar pentagon symbol letter in [1807.09812]*)


(* ::Text::Initialization:: *)
(*The two-loop pentagon symbol letter provided in Gehrmann, Henn and Presti's paper: "Pentagon functions for massless planar scattering amplitudes".*)


Wlist={v[1],v[2],v[3],v[4],v[5],v[3]+v[4],v[4]+v[5],v[1]+v[5],
v[1]+v[2],v[2]+v[3],v[1]-v[4],v[2]-v[5],-v[1]+v[3],-v[2]+v[4],
-v[3]+v[5],-v[1]-v[2]+v[4],-v[2]-v[3]+v[5],v[1]-v[3]-v[4],
v[2]-v[4]-v[5],-v[1]+v[3]-v[5],-v[1]-v[2]+v[3]+v[4],-v[2]-v[3]+v[4]+v[5],
v[1]-v[3]-v[4]+v[5],v[1]+v[2]-v[4]-v[5],-v[1]+v[2]+v[3]-v[5],
(-Sqrt[d]+v[1] v[2]-v[2] v[3]+v[3] v[4]-v[1] v[5]-v[4] v[5])/(Sqrt[d]+v[1] v[2]-v[2] v[3]+v[3] v[4]-v[1] v[5]-v[4] v[5]),
(-Sqrt[d]-v[1] v[2]+v[2] v[3]-v[3] v[4]-v[1] v[5]+v[4] v[5])/(Sqrt[d]-v[1] v[2]+v[2] v[3]-v[3] v[4]-v[1] v[5]+v[4] v[5]),
(-Sqrt[d]-v[1] v[2]-v[2] v[3]+v[3] v[4]+v[1] v[5]-v[4] v[5])/(Sqrt[d]-v[1] v[2]-v[2] v[3]+v[3] v[4]+v[1] v[5]-v[4] v[5]),
(-Sqrt[d]+v[1] v[2]-v[2] v[3]-v[3] v[4]-v[1] v[5]+v[4] v[5])/(Sqrt[d]+v[1] v[2]-v[2] v[3]-v[3] v[4]-v[1] v[5]+v[4] v[5]),
(-Sqrt[d]-v[1] v[2]+v[2] v[3]-v[3] v[4]+v[1] v[5]-v[4] v[5])/(Sqrt[d]-v[1] v[2]+v[2] v[3]-v[3] v[4]+v[1] v[5]-v[4] v[5])
,Sqrt[d]}/.d:>(v[1]v[2]+v[2]v[3]-v[3]v[4]+v[4]v[5]-v[5]v[1])^2-4 v[1]v[2]v[3](v[2]-v[4]-v[5]);


(* ::Section::Initialization::Closed:: *)
(*Evaluate spinor helicity bracket and Plucker coordinate*)


(* ::Subsection::Initialization::Closed:: *)
(*Definitions*)


(* ::Input::Initialization:: *)
ZZ=Table[x[i,j],{i,1,5},{j,1,3}];


(* ::Input::Initialization:: *)
randEval={x[1,1]->167,x[1,2]->241,x[1,3]->233,x[2,1]->433,x[2,2]->109,x[2,3]->79,x[3,1]->151,x[3,2]->173,x[3,3]->359,x[4,1]->97,x[4,2]->101,x[4,3]->53,x[5,1]->67,x[5,2]->467,x[5,3]->401};


(* ::Input::Initialization:: *)
(*randEval=Table[x[i,j]->RandomPrime[{50,500}],{i,1,6},{j,1,4}]//Flatten*)


(* ::Input::Initialization:: *)
sToP:={SpinorAngle[a_,b_]:>Signature[{a,b}]P@@Sort[{a,b}],SpinorSquare[a_,b_]:>Signature[{a,b}](-1)^(a+b+1) P@@Complement[{1,2,3,4,5},{a,b}]}


(* ::Input::Initialization:: *)
evalP:={P[a__]:>(B[a]/.randEval)}


(* ::Input::Initialization:: *)
B[a_,b_,c_]:=Det[{ZZ[[a]],ZZ[[b]],ZZ[[c]]}]
B[a_,b_]:=Det[{ZZ[[a,1;;2]],ZZ[[b,1;;2]]}]


(* ::Input::Initialization:: *)
penkin:={v[1]->SpinorAngle[1,2] SpinorSquare[1,2],v[2]->SpinorAngle[2,3] SpinorSquare[2,3],v[3]->SpinorAngle[3,4] SpinorSquare[3,4],v[4]->SpinorAngle[4,5] SpinorSquare[4,5],v[5]->SpinorAngle[5,1] SpinorSquare[5,1]}


(* ::Subsection::Initialization::Closed:: *)
(*Examples*)


(* ::Text::Initialization:: *)
(*hexkin maps the variables in the planar data (v[i],r[i],eps[i,j,k,l]) to spinor helicity bracket*)


(* ::Input::Initialization:: *)
v[1]/.penkin


(* ::Text::Initialization:: *)
(*sToP maps spinor helicity bracket to Plucker coordinate of Subscript[Fl, 2,4;6]*)


(* ::Input::Initialization:: *)
v[1]/.penkin/.sToP


(* ::Text::Initialization:: *)
(*evalP evaluate the numerical value of the Plucker coordinate with respect to a fixed random matrix*)


(* ::Input::Initialization:: *)
v[1]/.penkin/.sToP/.evalP


(* ::Subsection::Initialization::Closed:: *)
(*Sanity checks: antisymmetry of spinor bracket, Shouten identity and momentum conservation*)


(* ::Input::Initialization:: *)
Table[SpinorAngle[i,j]+SpinorAngle[j,i],{i,1,5},{j,1,5}]/.sToP/.evalP//Flatten//DeleteDuplicates


(* ::Input::Initialization:: *)
Table[SpinorSquare[i,j]+SpinorSquare[j,i],{i,1,5},{j,1,5}]/.sToP/.evalP//Flatten//DeleteDuplicates


(* ::Input::Initialization:: *)
Table[SpinorAngle[i,j]SpinorAngle[k,l]+SpinorAngle[j,k]SpinorAngle[i,l]+SpinorAngle[k,i]SpinorAngle[j,l],{i,1,5},{j,1,5},{k,1,5},{l,1,5}]/.sToP/.evalP//Flatten//DeleteDuplicates


(* ::Input::Initialization:: *)
Table[Sum[SpinorAngle[i,j]SpinorSquare[j,k],{j,1,5}],{i,1,5},{k,1,5}]/.sToP/.evalP//Flatten//DeleteDuplicates


(* ::Section::Initialization::Closed:: *)
(*Cluster variables of Subscript[Fl, 2,3;5]*)


cyclic5=Table[Thread[{1,2,3,4,5}->i],{i,{{1,2,3,4,5},{2,3,4,5,1},{3,4,5,1,2},{4,5,1,2,3},{5,1,2,3,4}}}];


{a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10]}=
(P[1,2]/.cyclic5)~Join~(P[1,3]/.cyclic5);

{b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10]}=
(P[3,4,5]/.cyclic5)~Join~(P[2,4,5]/.cyclic5);

{c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9],c[10],c[11],c[12],c[13],c[14],c[15]}=
(-P[3,5] P[1,2,4]-P[1,2] P[3,4,5]/.cyclic5)~Join~
(P[4,5] P[1,2,3]-P[1,2] P[3,4,5]/.cyclic5)~Join~
(P[2,5] P[1,3,4]-P[1,4] P[2,3,5]/.cyclic5);


(* ::Section::Initialization::Closed:: *)
(*Matching symbol letters*)


(* ::Text::Initialization:: *)
(*We provide the combination of the cluster variables to match with the letters.*)
(*Note that the replacement rules are only valid for Re[Subscript[\[Epsilon], 1,2,3,4]]>0. *)
(*For Re[Subscript[\[Epsilon], 1,2,3,4]]<0, Subscript[W, 26]-Subscript[W, 30] map to the inverse of what is provided below.*)


(*Subscript[\[Epsilon], 1,2,3,4]*)
SpinorAngle[1,2]SpinorAngle[4,5]SpinorSquare[2,4]SpinorSquare[1,5]-
SpinorSquare[1,2]SpinorSquare[4,5]SpinorAngle[2,4]SpinorAngle[1,5]/.sToP/.evalP


{W[1]->a[1] b[1],W[2]->a[2] b[2],W[3]->a[3] b[3],W[4]->a[4] b[4],W[5]->a[5] b[5],
W[6]->-c[1],W[7]->-c[2],W[8]->-c[3],W[9]->-c[4],W[10]->-c[5],W[11]->-c[6],W[12]->-c[7],
W[13]->-c[8],W[14]->-c[9],W[15]->-c[10],W[16]->-a[6] b[6],W[17]->-a[7] b[7],
W[18]->-a[8] b[8],W[19]->-a[9] b[9],W[20]->-a[10] b[10],W[21]->c[11],W[22]->c[12],
W[23]->c[13],W[24]->c[14],W[25]->c[15],W[26]->-((a[1] a[4] b[5] b[7])/(a[5] a[7] b[1] b[4])),
W[27]->-((a[2] a[5] b[1] b[8])/(a[1] a[8] b[2] b[5])),
W[28]->-((a[1] a[3] b[2] b[9])/(a[2] a[9] b[1] b[3])),
W[29]->-((a[2] a[4] b[3] b[10])/(a[3] a[10] b[2] b[4])),
W[30]->-((a[3] a[5] b[4] b[6])/(a[4] a[6] b[3] b[5]))}


Table[Wlist[[i]],{i,1,30}]/.penkin/.sToP/.evalP;
{a[1]b[1],a[2]b[2],a[3]b[3],a[4]b[4],a[5]b[5],-c[1],-c[2],-c[3],-c[4],-c[5],
-c[6],-c[7],-c[8],-c[9],-c[10],-a[6]b[6],-a[7]b[7],-a[8]b[8],-a[9]b[9],-a[10]b[10],
c[11],c[12],c[13],c[14],c[15],-((a[1]a[4]b[7]b[5])/(a[7]a[5]b[1]b[4])),-((a[2]a[5]b[8]b[1])/(a[8]a[1]b[2]b[5])),
-((a[3]a[1]b[9]b[2])/(a[9]a[2]b[3]b[1])),-((a[4]a[2]b[10]b[3])/(a[10]a[3]b[4]b[2])),-((a[5]a[3]b[6]b[4])/(a[6]a[4]b[5]b[3]))
}/.evalP;
%-%%
