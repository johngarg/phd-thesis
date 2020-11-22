(* ::Package:: *)

Subscript[PMNS\[Theta], 12]=33.36 Degree;
Subscript[PMNS\[Theta], 23]=40.4 Degree;
Subscript[PMNS\[Theta]IO, 23]=50 Degree;
Subscript[PMNS\[Theta], 13]=8.66 Degree;
Subscript[PMNS\[Delta], CP]=0 Degree;
PMNS=({
 {Cos[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 13]], Sin[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 13]], Sin[Subscript[PMNS\[Theta], 13]]E^(-I Subscript[PMNS\[Delta], CP])},
 {-Sin[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 23]]-Cos[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta], 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), Cos[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 23]]-Sin[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta], 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), Sin[Subscript[PMNS\[Theta], 23]]Cos[Subscript[PMNS\[Theta], 13]]},
 {Sin[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta], 23]]-Cos[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), -Cos[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta], 23]]-Sin[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), Cos[Subscript[PMNS\[Theta], 23]]Cos[Subscript[PMNS\[Theta], 13]]}
});
PMNSIO=({
 {Cos[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 13]], Sin[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta], 13]], Sin[Subscript[PMNS\[Theta], 13]]E^(-I Subscript[PMNS\[Delta], CP])},
 {-Sin[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta]IO, 23]]-Cos[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta]IO, 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), Cos[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta]IO, 23]]-Sin[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta]IO, 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), Sin[Subscript[PMNS\[Theta]IO, 23]]Cos[Subscript[PMNS\[Theta], 13]]},
 {Sin[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta]IO, 23]]-Cos[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta]IO, 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), -Cos[Subscript[PMNS\[Theta], 12]]Sin[Subscript[PMNS\[Theta]IO, 23]]-Sin[Subscript[PMNS\[Theta], 12]]Cos[Subscript[PMNS\[Theta]IO, 23]]Sin[Subscript[PMNS\[Theta], 13]]E^(I Subscript[PMNS\[Delta], CP]), Cos[Subscript[PMNS\[Theta]IO, 23]]Cos[Subscript[PMNS\[Theta], 13]]}
});


massEigen=10^-9 {0,Sqrt[7.5 10^-5],Sqrt[2.473 10^-3]};
massEigenIO=10^-9 {Sqrt[2.514 10^-3-7.5 10^-5],Sqrt[2.514 10^-3],0};  (* \[Nu] mass in GeV, index: lepton *)

ghat[s_,t_]:=Block[{w,xplus,xminus},
w=Sqrt[1+s^2+t^2-2(s+t+s t)];xplus=1/2 (-1+s+t+w);xminus=1/2 (-1+s+t-w);s/2 Log[s]Log[t] + (s(1-s)+3s t + 2(1-t)xplus)/(2w) (PolyLog[2,xplus/(xplus-s)]-PolyLog[2,(xplus-s)/xplus]+PolyLog[2,(t-1)/xplus]-PolyLog[2,(t-1)/(xplus-s)])-(s(1-s)+3s t + 2(1-t)xminus)/(2w) (PolyLog[2,xminus/(xminus-s)]-PolyLog[2,(xminus-s)/xminus]+PolyLog[2,(t-1)/xminus]-PolyLog[2,(t-1)/(xminus-s)])];
ghat[s_,0]:=-s \[Pi]^2/6 - (1-s)Log[s]Log[1-s]-(1-s)PolyLog[2,s];
(* Integral under approximation that new particles much heavier than down type quarks *)
Subscript[Int, \[Alpha]_,\[Beta]_][mf_,m\[Phi]_]:=Module[{t},
t[n_]:=m\[Phi][[n]]^2/mf^2;
\[Pi]^4/mf^2 (ghat[t[\[Alpha]],0]-ghat[t[\[Alpha]],t[\[Beta]]])/(t[\[Alpha]]t[\[Beta]])];

IntMatrix[mf_,m\[Phi]_]:=Array[Subscript[Int, #1,#2][mf,m\[Phi]]&,{2,2}];


RMatrix[\[Theta]_]:=({
 {0, Cos[\[Theta]], Sin[\[Theta]]},
 {0, -Sin[\[Theta]], Cos[\[Theta]]}
})//Transpose; (* indices: lepton, LQ *)
RMatrixIO[\[Theta]_]:=({
 {Cos[\[Theta]], Sin[\[Theta]], 0},
 {-Sin[\[Theta]], Cos[\[Theta]], 0}
})//Transpose; (* indices: lepton, LQ *)


IntInvRoot[mf_,m\[Phi]_]:=MatrixPower[IntMatrix[mf,m\[Phi]]//N,-0.5];

\[CapitalOmega]Matrix[mf_,m\[Phi]_,\[Theta]_]:= MatrixPower[DiagonalMatrix[massEigen],0.5].RMatrix[\[Theta]].IntInvRoot[mf,m\[Phi]];
\[CapitalOmega]MatrixIO[mf_,m\[Phi]_,\[Theta]_]:= MatrixPower[DiagonalMatrix[massEigenIO],0.5].RMatrixIO[\[Theta]].IntInvRoot[mf,m\[Phi]];
\[CapitalLambda][mf_,m\[Phi]_,\[Theta]_]:=With[{mb=4.2},(2\[Pi])^4/(2mb Sqrt[mf]) PMNS.\[CapitalOmega]Matrix[mf,m\[Phi],\[Theta]] (* indices: lepton, LQ *)];
\[CapitalLambda]IO[mf_,m\[Phi]_,\[Theta]_]:=With[{mb=4.2},(2\[Pi])^4/(2mb Sqrt[mf]) PMNSIO.\[CapitalOmega]MatrixIO[mf,m\[Phi],\[Theta]] (* indices: lepton, LQ *)];
minimalm\[Nu][mf_,m\[Phi]_,\[Theta]_]:=With[{mb=4.2},4 (mf mb^2)/(2\[Pi])^8 \[CapitalLambda][mf,m\[Phi],\[Theta]].(IntMatrix[mf,m\[Phi]]).Transpose[\[CapitalLambda][mf,m\[Phi],\[Theta]]]];


Br\[Mu]NeN[m\[Phi]_,\[Lambda]LeuMatrix_,\[Lambda]ReuMatrix_]:=Module[{t,gZLVd,gZLVu,gPLV,gNLV,A2L,A2R,\[Sigma]L,\[Sigma]R,g\[Gamma]LVd,g\[Gamma]RVd,g\[Gamma]LVu,g\[Gamma]RVu,gZRVd,gZRVu,gLVu,gLVd,gRVu,gRVd,gPRV,gNRV,\[Omega]conv,\[Alpha]=1/137.},
Subscript[t, q_]:=mu[[q]]^2/m\[Phi]^2;
g\[Gamma]LVd=Sum[-(\[Alpha]/(144\[Pi])) ( \[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]])/m\[Phi]^2 ((Subscript[t, a]^3-18Subscript[t, a]^2+27Subscript[t, a]+2(Subscript[t, a]^3+6Subscript[t, a]-4)Log[Subscript[t, a]]-10)/(Subscript[t, a]-1)^4),{a,1,3}];
g\[Gamma]RVd=Sum[-(\[Alpha]/(144\[Pi])) ( \[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]])/m\[Phi]^2 ((Subscript[t, a]^3-18Subscript[t, a]^2+27Subscript[t, a]+2(Subscript[t, a]^3+6Subscript[t, a]-4)Log[Subscript[t, a]]-10)/(Subscript[t, a]-1)^4),{a,1,3}];
gZLVd=Sum[-(8/Sqrt[2])GF (4Sin[\[Theta]w]^2-3)/(128 \[Pi]^2) \[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]]((Subscript[t, a](Subscript[t, a]-Log[Subscript[t, a]]-1))/(Subscript[t, a]-1)^2),{a,1,3}];
gZRVd=Sum[-(8/Sqrt[2])GF (4Sin[\[Theta]w]^2-3)/(128 \[Pi]^2) \[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]]((Subscript[t, a](Subscript[t, a]-Log[Subscript[t, a]]-1))/(Subscript[t, a]-1)^2),{a,1,3}];
g\[Gamma]LVu=-2g\[Gamma]LVd;
g\[Gamma]RVu=-2g\[Gamma]RVd;
gZRVu=-((8Sin[\[Theta]w]^2-3)/(4Sin[\[Theta]w]^2-3))gZRVd;
gZLVu=-((8Sin[\[Theta]w]^2-3)/(4Sin[\[Theta]w]^2-3))gZLVd;
gLVu=gZLVu+g\[Gamma]LVu;
gLVd=gZLVd+g\[Gamma]LVd;
gRVu=gZRVu+g\[Gamma]RVu;
gRVd=gZRVd+g\[Gamma]RVd;

gPLV=2gLVu+gLVd;
gNLV=gLVu+2gLVd;

gPRV=2gRVu+gRVd;
gNRV=gRVu+2gRVd;

\[Sigma]L=Sum[1/(16\[Pi]^2) (\[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]]m\[Mu])/m\[Phi]^2 (1+4Subscript[t, a]-5Subscript[t, a]^2+2Subscript[t, a](2+Subscript[t, a])Log[Subscript[t, a]])/(4(Subscript[t, a]-1)^4)+1/(16\[Pi]^2) (\[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]]mu[[a]])/m\[Phi]^2 (7-8Subscript[t, a]+Subscript[t, a]^2+2(2+Subscript[t, a])Log[Subscript[t, a]])/(2(Subscript[t, a]-1)^3),{a,1,3}];
\[Sigma]R=Sum[1/(16\[Pi]^2) (\[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]]m\[Mu])/m\[Phi]^2 (1+4Subscript[t, a]-5Subscript[t, a]^2+2Subscript[t, a](2+Subscript[t, a])Log[Subscript[t, a]])/(4(Subscript[t, a]-1)^4)+1/(16\[Pi]^2) (\[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]]mu[[a]])/m\[Phi]^2 (7-8Subscript[t, a]+Subscript[t, a]^2+2(2+Subscript[t, a])Log[Subscript[t, a]])/(2(Subscript[t, a]-1)^3),{a,1,3}];
A2L=\[Sigma]L/m\[Mu];
A2R=\[Sigma]R/m\[Mu];

\[Omega]conv=4m\[Mu]^5 Abs[1/8 Conjugate[A2R]0.108 + gPLV 0.061 + gNLV 0.0859]^2+4m\[Mu]^5 Abs[1/8 Conjugate[A2L]0.108 + gPRV 0.061 + gNRV 0.0859]^2;
\[Omega]conv/(13.07 10^6 (4.135668*10^-24))
];
Br\[Mu]NeNCOMET[m\[Phi]_,\[Lambda]LeuMatrix_,\[Lambda]ReuMatrix_]:=Module[{t,gZLVd,gZLVu,gPLV,gNLV,A2L,A2R,\[Sigma]L,\[Sigma]R,g\[Gamma]LVd,g\[Gamma]RVd,g\[Gamma]LVu,g\[Gamma]RVu,gZRVd,gZRVu,gLVu,gLVd,gRVu,gRVd,gPRV,gNRV,\[Omega]conv,\[Alpha]=1/137.},
Subscript[t, q_]:=mu[[q]]^2/m\[Phi]^2;
g\[Gamma]LVd=Sum[-(\[Alpha]/(144\[Pi])) ( \[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]])/m\[Phi]^2 ((Subscript[t, a]^3-18Subscript[t, a]^2+27Subscript[t, a]+2(Subscript[t, a]^3+6Subscript[t, a]-4)Log[Subscript[t, a]]-10)/(Subscript[t, a]-1)^4),{a,1,3}];
g\[Gamma]RVd=Sum[-(\[Alpha]/(144\[Pi])) ( \[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]])/m\[Phi]^2 ((Subscript[t, a]^3-18Subscript[t, a]^2+27Subscript[t, a]+2(Subscript[t, a]^3+6Subscript[t, a]-4)Log[Subscript[t, a]]-10)/(Subscript[t, a]-1)^4),{a,1,3}];
gZLVd=Sum[-(8/Sqrt[2])GF (4Sin[\[Theta]w]^2-3)/(128 \[Pi]^2) \[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]]((Subscript[t, a](Subscript[t, a]-Log[Subscript[t, a]]-1))/(Subscript[t, a]-1)^2),{a,1,3}];
gZRVd=Sum[-(8/Sqrt[2])GF (4Sin[\[Theta]w]^2-3)/(128 \[Pi]^2) \[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]]((Subscript[t, a](Subscript[t, a]-Log[Subscript[t, a]]-1))/(Subscript[t, a]-1)^2),{a,1,3}];
g\[Gamma]LVu=-2g\[Gamma]LVd;
g\[Gamma]RVu=-2g\[Gamma]RVd;
gZRVu=-((8Sin[\[Theta]w]^2-3)/(4Sin[\[Theta]w]^2-3))gZRVd;
gZLVu=-((8Sin[\[Theta]w]^2-3)/(4Sin[\[Theta]w]^2-3))gZLVd;
gLVu=gZLVu+g\[Gamma]LVu;
gLVd=gZLVd+g\[Gamma]LVd;
gRVu=gZRVu+g\[Gamma]RVu;
gRVd=gZRVd+g\[Gamma]RVd;

gPLV=2gLVu+gLVd;
gNLV=gLVu+2gLVd;

gPRV=2gRVu+gRVd;
gNRV=gRVu+2gRVd;

\[Sigma]L=Sum[1/(16\[Pi]^2) (\[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]]m\[Mu])/m\[Phi]^2 (1+4Subscript[t, a]-5Subscript[t, a]^2+2Subscript[t, a](2+Subscript[t, a])Log[Subscript[t, a]])/(4(Subscript[t, a]-1)^4)+1/(16\[Pi]^2) (\[Lambda]LeuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]]mu[[a]])/m\[Phi]^2 (7-8Subscript[t, a]+Subscript[t, a]^2+2(2+Subscript[t, a])Log[Subscript[t, a]])/(2(Subscript[t, a]-1)^3),{a,1,3}];
\[Sigma]R=Sum[1/(16\[Pi]^2) (\[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]ReuMatrix[[1,a]]]m\[Mu])/m\[Phi]^2 (1+4Subscript[t, a]-5Subscript[t, a]^2+2Subscript[t, a](2+Subscript[t, a])Log[Subscript[t, a]])/(4(Subscript[t, a]-1)^4)+1/(16\[Pi]^2) (\[Lambda]ReuMatrix[[2,a]]Conjugate[\[Lambda]LeuMatrix[[1,a]]]mu[[a]])/m\[Phi]^2 (7-8Subscript[t, a]+Subscript[t, a]^2+2(2+Subscript[t, a])Log[Subscript[t, a]])/(2(Subscript[t, a]-1)^3),{a,1,3}];
A2L=\[Sigma]L/m\[Mu];
A2R=\[Sigma]R/m\[Mu];

(*Print[{"g\[Gamma]LVd"->g\[Gamma]LVd,"g\[Gamma]RVd"->g\[Gamma]RVd,"gZLVd"->gZLVd,"gZRVd"->gZRVd,"gLVu"\[Rule]gLVu,"gLVd"->gLVd,"gRVu"->gRVu,"gRVd"->gRVd,"A2L"->A2L,"A2R"->A2R}];*)

\[Omega]conv=4m\[Mu]^5 Abs[1/8 Conjugate[A2R]0.0159 + gPLV 0.0163 + gNLV 0.0357]^2+4m\[Mu]^5 Abs[1/8 Conjugate[A2L]0.0159 + gPRV 0.0163 + gNRV 0.0357]^2;
\[Omega]conv/(0.7054 10^6 (4.135668*10^-24))
];


\[Mu]NeN[m_,\[Theta]_,w3_]:=Br\[Mu]NeN[m,({
 {0, 0, Chop[Transpose[\[CapitalLambda][25000,{m,20000},\[Theta]]][[1]]][[1]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda][25000,{m,20000},\[Theta]]][[1]]][[2]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda][25000,{m,20000},\[Theta]]][[1]]][[3]]/w3}
}).ConjugateTranspose[CKM],({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];
\[Mu]NeNCOMET[m_,\[Theta]_,w3_]:=Br\[Mu]NeNCOMET[m,({
 {0, 0, Chop[Transpose[\[CapitalLambda][25000,{m,10000},\[Theta]]][[1]]][[1]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda][25000,{m,10000},\[Theta]]][[1]]][[2]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda][25000,{m,10000},\[Theta]]][[1]]][[3]]/w3}
}).ConjugateTranspose[CKM],({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];

\[Mu]NeNIO[m_,\[Theta]_,w3_]:=Br\[Mu]NeN[m,({
 {0, 0, Chop[Transpose[\[CapitalLambda]IO[25000,{m,20000},\[Theta]]][[1]]][[1]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda]IO[25000,{m,20000},\[Theta]]][[1]]][[2]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda]IO[25000,{m,20000},\[Theta]]][[1]]][[3]]/w3}
}).ConjugateTranspose[CKM],({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];
\[Mu]NeNCOMETIO[m_,\[Theta]_,w3_]:=Br\[Mu]NeNCOMET[m,({
 {0, 0, Chop[Transpose[\[CapitalLambda]IO[25000,{m,10000},\[Theta]]][[1]]][[1]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda]IO[25000,{m,10000},\[Theta]]][[1]]][[2]]/w3},
 {0, 0, Chop[Transpose[\[CapitalLambda]IO[25000,{m,10000},\[Theta]]][[1]]][[3]]/w3}
}).ConjugateTranspose[CKM],({
 {0, 0, 0},
 {0, 0, 0},
 {0, 0, 0}
})];
