(* ::Package:: *)

Print["FlavorConstraintsLQ Version: 1.0 (8 Feb 2017)."];
Print["Author: John Gargalionis"];
Print["Please cite arXiv:1702.xxxxx and the appropriate references therein."];
Print["Note that the package FlavorConstraintsLQ needs to be loaded after PhysicsWorkbench and BAnomalies."];


BeginPackage["FlavorConstraintsLQ`"];

{(* What is this section actually for? *)}

EndPackage[];


<<"~/Dropbox/mathematicaPackages/RunDec.m"

\[Alpha]sMt5 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mt/.NumDef,5,3];
\[Alpha]sMt6 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mt/.NumDef,6,3];
\[Alpha]sMb = AlphasExact[asMz/.NumDef,Mz/.NumDef,4.2,5,3];
\[Alpha]sMW = AlphasExact[asMz/.NumDef,Mz/.NumDef,80.4,5,3];
\[Alpha]sMLQ[m_] := AsRunDec[asMz/.NumDef,Mz/.NumDef,m,3];

(* running factors for scalar and tensor operators *)
runningS[m_]:=(\[Alpha]sMt5/\[Alpha]sMb)^(-12/23)(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(-12/21);
runningT[m_]:=(\[Alpha]sMt5/\[Alpha]sMb)^(4/23)(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(4/21);

(* running for semi-leptonic scalar operators *)
\[Alpha]sMb4 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mb/.NumDef,4,3];
\[Alpha]sMb5 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mb/.NumDef,5,3];
\[Alpha]sMc3 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mc/.NumDef,3,3];
\[Alpha]sMc4 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mc/.NumDef,4,2];
\[Alpha]s2 = AlphasExact[asMz/.NumDef,Mz/.NumDef,2.,4,2];
beta0[nf_]:=11-2 nf/3;
runningSK[m_]:=With[{\[Gamma]=-8.},(\[Alpha]sMb4/\[Alpha]s2)^(\[Gamma]/(2 beta0[4]))(\[Alpha]sMt5/\[Alpha]sMb5)^(\[Gamma]/(2 beta0[5]))(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(\[Gamma]/(2 beta0[6]))];
runningSD[m_]:=With[{\[Gamma]=-8.},(\[Alpha]sMb4/\[Alpha]sMc4)^(\[Gamma]/(2 beta0[4]))(\[Alpha]sMt5/\[Alpha]sMb5)^(\[Gamma]/(2 beta0[5]))(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(\[Gamma]/(2 beta0[6]))];
runningTD[m_]:=With[{\[Gamma]=8./3},(\[Alpha]sMb4/\[Alpha]sMc4)^(\[Gamma]/(2 beta0[4]))(\[Alpha]sMt5/\[Alpha]sMb5)^(\[Gamma]/(2 beta0[5]))(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(\[Gamma]/(2 beta0[6]))];


(* Bauer & Neubert *)
rK\[Nu]\[Nu]Condition[m_,x_]:=With[{rr=1.91/(m/1000)^2},
(1-(2 rr)/3 Re[
Sum[Conjugate[x[[i,2]]]x[[i,3]],{i,3}]/(CKM[[3,3]]CKM[[3,2]])
]+ rr^2/3 (Sum[x[[i,3]]Conjugate[x[[i,3]]],{i,3}]Sum[x[[i,2]]Conjugate[x[[i,2]]],{i,3}])/Abs[CKM[[3,3]]Conjugate[CKM[[3,2]]]]^2)<4.3];

cBsCondition[m_,x_]:=Block[{\[Alpha]=1/127.,running},
running = (\[Alpha]sMt5/\[Alpha]sMW)^(6/23)(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(2/7);
Abs[1+running 1/((Sqrt[4\[Pi] \[Alpha]]/Sqrt[0.23129])^4 2.30) (mW^2/(m)^2)(Sum[Conjugate[x[[i,2]]]x[[i,3]],{i,3}]/(CKM[[3,3]]CKM[[3,2]]))^2]<(1.052+2 0.084)];

brD\[Mu]\[Mu][m_,x_,y_]:=Block[{fD,\[Tau]D=1040 10^-15,\[CapitalGamma]D,mc=1.27,\[Beta],z},
fD=212 10^-3;
\[CapitalGamma]D=hbar/\[Tau]D;
z=x.ConjugateTranspose[CKM];
\[Beta]=Sqrt[1-4 m\[Mu]^2/mD^2];
1/\[CapitalGamma]D (fD^2 mD^3)/(512 \[Pi] m^4) (mD/mc)^2 \[Beta] (
\[Beta]^2 Abs[runningSD[m]]^2 Abs[-z[[2,2]] y[[2,1]] +Conjugate[ y[[2,2]]] Conjugate[z[[2,1]]]]^2
+ Abs[
runningSD[m](-z[[2,2]]y[[2,1]] -Conjugate[y[[2,2]]] Conjugate[z[[2,1]]])+(2m\[Mu] mc)/mD^2 (z[[2,2]]Conjugate[z[[2,1]]]+Conjugate[y[[2,2]]]y[[2,1]])]^2)
];
brD\[Mu]\[Mu]Condition[m_,x_,y_]:=brD\[Mu]\[Mu][m,x,y]<7.6 10^-9;


(* Becirevic et al. *)

brB\[Tau]\[Nu]Condition[m_,x_,y_]:=With[{fB=186 10^-3,\[Tau]B=1.638 10^-12,z=x.ConjugateTranspose[CKM],mu=0.0022,mb=4.18},
(1.14-2 0.27)10^-4<\[Tau]B/hbar ((GF^2 mB Abs[CKM[[1,3]]]^2)/(8\[Pi]) fB^2 m\[Tau]^2 (1-m\[Tau]^2/mB^2)^2 Sum[Abs[KroneckerDelta[3,i]+1/(4Sqrt[2]GF CKM[[1,3]]) (Conjugate[z[[3,1]]]x[[i,3]])/m^2-1/(4Sqrt[2]GF CKM[[1,3]]) (runningS[m]y[[3,1]]x[[i,3]])/m^2 (mB^2/(m\[Tau](mu+mb)))]^2,{i,3}])<(1.14+2 0.27)10^-4
];

brB\[Mu]\[Nu]Condition[m_,x_,y_]:=With[{fB=186 10^-3,\[Tau]B=1.638 10^-12,z=x.ConjugateTranspose[CKM],mu=0.0022,mb=4.18},
\[Tau]B/hbar ((GF^2 mB Abs[CKM[[1,3]]]^2)/(8\[Pi]) fB^2 mLeps[[2]]^2 (1-mLeps[[2]]^2/mB^2)^2 Sum[Abs[KroneckerDelta[2,i]+1/(4Sqrt[2]GF CKM[[1,3]]) (Conjugate[z[[2,1]]]x[[i,3]])/m^2-1/(4Sqrt[2]GF CKM[[1,3]]) (runningS[m]y[[2,1]]x[[i,3]])/m^2 (mB^2/(m\[Mu](mu+mb)))]^2,{i,3}])<10^-6
];

brK\[Mu]\[Nu]Condition[m_,x_,y_]:=With[{fK=155.6 10^-3,\[Tau]K=1.238 10^-8,z=x.ConjugateTranspose[CKM],mK=493.677 10^-3,mu=0.0022},
(63.56-2 0.11)10^-2<\[Tau]K/hbar((GF^2mK Abs[CKM[[1,2]]]^2)/(8\[Pi])fK^2m\[Mu]^2(1-m\[Mu]^2/mK^2)^2Sum[Abs[KroneckerDelta[2,i]+1/(4Sqrt[2]GF CKM[[1,2]])(Conjugate[z[[2,1]]]x[[i,2]])/m^2-1/(4Sqrt[2]GF CKM[[1,2]])(runningSK[m]y[[2,1]]x[[i,2]])/m^2(mK^2/(m\[Mu](mu+md[[2]])))]^2,{i,3}])<(63.56+2 0.11)10^-2
];

brDs\[Mu]\[Nu]Condition[m_,x_,y_]:=With[{fDs=249 10^-3,\[Tau]Ds=500 10^-15,z=x.ConjugateTranspose[CKM],mDs=1.9683,mc=1.27,ms=0.095},
(5.56-2 0.25)10^-3<\[Tau]Ds/hbar ((GF^2 mDs Abs[CKM[[2,2]]]^2)/(8\[Pi]) fDs^2 m\[Mu]^2 (1-m\[Mu]^2/mDs^2)^2 Sum[Abs[KroneckerDelta[2,i]+1/(4Sqrt[2]GF CKM[[2,2]]) (Conjugate[z[[2,2]]]x[[i,2]])/m^2-1/(4Sqrt[2]GF CKM[[2,2]]) (runningSD[m]y[[2,2]]x[[i,2]])/m^2 (mDs^2/(m\[Mu](mc+ms)))]^2,{i,3}])<(5.56+2 0.25)10^-3
];

brDs\[Tau]\[Nu]Condition[m_,x_,y_]:=With[{fDs=249 10^-3,\[Tau]Ds=500 10^-15,z=x.ConjugateTranspose[CKM],mDs=1.9683,mc=1.27,ms=0.095},
(5.55-2 0.24)10^-2<\[Tau]Ds/hbar ((GF^2 mDs Abs[CKM[[2,2]]]^2)/(8\[Pi]) fDs^2 m\[Tau]^2 (1-m\[Tau]^2/mDs^2)^2 Sum[Abs[KroneckerDelta[3,i]+1/(4Sqrt[2]GF CKM[[2,2]]) (Conjugate[z[[3,2]]]x[[i,2]])/m^2-1/(4Sqrt[2]GF CKM[[2,2]]) (runningSD[m]y[[3,2]]x[[i,2]])/m^2 (mDs^2/(m\[Tau](mc+ms)))]^2,{i,3}])<(5.55+2 0.24)10^-2
];

(* new processes and ratios *)

rKe\[Mu]Condition[m_,x_,y_]:=Block[{fK=155.6 10^-3,\[Tau]K=1.238 10^-8,z=x.ConjugateTranspose[CKM],mK=493.677 10^-3,mu=0.0022,\[CapitalGamma]K\[Mu]\[Nu],\[CapitalGamma]Ke\[Nu]},\[CapitalGamma]K\[Mu]\[Nu]=(m\[Mu]^2(1-m\[Mu]^2/mK^2)^2Sum[Abs[KroneckerDelta[2,i]+1/(4Sqrt[2]GF CKM[[1,2]])(Conjugate[z[[2,1]]]x[[i,2]])/m^2-1/(4Sqrt[2]GF CKM[[1,2]])(runningSK[m]y[[2,1]]x[[i,2]])/m^2(mK^2/(m\[Mu](mu+md[[2]])))]^2,{i,3}]);
\[CapitalGamma]Ke\[Nu]=(me^2(1-me^2/mK^2)^2Sum[Abs[KroneckerDelta[1,i]]^2,{i,3}]);
(2.488-2 0.009) 10^-5<((\[CapitalGamma]Ke\[Nu]/\[CapitalGamma]K\[Mu]\[Nu] )(1- 0.03786 + 0.00135)(1+0.0055))<(2.488+2 0.009) 10^-5
];

rK\[Tau]\[Mu]Condition[m_,x_,y_]:=Module[{z,fK=155.6 10^-3,\[Tau]\[Tau]=2.9 10^-13,gL,gRR,\[Tau]K=1.238 10^-8},
z=x.ConjugateTranspose[CKM];
gL[k_]:=vev^2/(2m^2) z[[3,1]] x[[k,2]];
gRR[k_]:=-(vev^2/(2m^2))runningSK[m]x[[k,2]]y[[3,1]];

(1.101-2 0.016) 10^-2<(hbar/\[Tau]\[Tau])^-1 m\[Tau]/(16\[Pi]) (1-mK^2/m\[Tau]^2)^2 GF^2 fK^2 mK^2 Sum[Abs[m\[Tau]/mK (CKM[[1,2]] KroneckerDelta[3,k]+gL[k])+mK/(mu[[1]]+md[[2]]) gRR[k]]^2,{k,1,3}](\[Tau]K/hbar((GF^2mK Abs[CKM[[1,2]]]^2)/(8\[Pi])fK^2m\[Mu]^2(1-m\[Mu]^2/mK^2)^2Sum[Abs[KroneckerDelta[2,i]+1/(4Sqrt[2]GF CKM[[1,2]])(Conjugate[z[[2,1]]]x[[i,2]])/m^2-1/(4Sqrt[2]GF CKM[[1,2]])(y[[2,1]]x[[i,2]])/m^2(mK^2/(m\[Mu](mu[[1]]+md[[2]])))]^2,{i,3}]))^-1<(1.101+2 0.016) 10^-2
];


(* Expression from Bauer & Neubert, EW fit from arXiv:xxx.xxxxx *)
\[Delta]gL[m\[Phi]_,c1_,c2_,c3_]:=With[{mt=160},3/(32\[Pi]^2) mt^2/m\[Phi]^2 (Log[m\[Phi]^2/mt^2]-1)c3^2-1/(32\[Pi]^2) mZ^2/m\[Phi]^2 (c1^2+c2^2)((1 - 4/3 (0.23142))(Log[m\[Phi]^2/mZ^2]+I \[Pi]+1/3)-(0.23142)/9)]
\[Delta]gR[m\[Phi]_,c1_,c2_,c3_]:=With[{mt=160},-(3/(32\[Pi]^2)) mt^2/m\[Phi]^2 (Log[m\[Phi]^2/mt^2]-1)c3^2-1/(32\[Pi]^2) mZ^2/m\[Phi]^2 (c1^2+c2^2)((- (4/3)(0.23142))(Log[m\[Phi]^2/mZ^2]+I \[Pi]+1/3)-(0.23142)/9)]
zllConditionsEWFit[m_,x_,y_]:=With[{z=x.ConjugateTranspose[CKM]},
(-0.00085<Re[\[Delta]gL[m,z[[3,1]],z[[3,2]],z[[3,3]]]]<0.0012)&&(-0.00085<Re[\[Delta]gL[m,z[[2,1]],z[[2,2]],z[[2,3]]]]<0.0012)&&
(-0.00054<Re[\[Delta]gR[m,y[[3,1]],y[[3,2]],y[[3,3]]]]<0.00067)&&
(-0.00054<Re[\[Delta]gR[m,y[[2,1]],y[[2,2]],y[[2,3]]]]<0.00067)
];

zllConditions[m_,x_,y_]:=Block[{
n=2,measured\[Mu]gLgR,measured\[Tau]gLgR,gASM=-0.50127,gVSM=-0.03712,gLSM,gRSM,z=x.ConjugateTranspose[CKM]
},

gLSM=0.5(gVSM-gASM);
gRSM=0.5(gVSM+gASM);

measured\[Mu]gLgR=With[{gV=-0.0367,\[Delta]gV=0.0023,\[Delta]gA=0.00054,gA=-0.50111},
{
0.5{(gV-gA-n Sqrt[\[Delta]gV^2+\[Delta]gA^2]),(gV-gA+n Sqrt[\[Delta]gV^2+\[Delta]gA^2])},0.5{(gV+gA-n Sqrt[\[Delta]gV^2+\[Delta]gA^2]),(gV+gA+n Sqrt[\[Delta]gV^2+\[Delta]gA^2])}
}];

measured\[Tau]gLgR=With[{gV=-0.0366,\[Delta]gV=0.001,\[Delta]gA=0.00064,gA=-0.50204},
{
0.5{(gV-gA-n Sqrt[\[Delta]gV^2+\[Delta]gA^2]),(gV-gA+n Sqrt[\[Delta]gV^2+\[Delta]gA^2])},0.5{(gV+gA-n Sqrt[\[Delta]gV^2+\[Delta]gA^2]),(gV+gA+n Sqrt[\[Delta]gV^2+\[Delta]gA^2])}
}];

measured\[Tau]gLgR[[1,1]]<gLSM+Re[\[Delta]gL[m,z[[3,1]],z[[3,2]],z[[3,3]]]]<measured\[Tau]gLgR[[1,2]]&&measured\[Mu]gLgR[[1,1]]<gLSM+Re[\[Delta]gL[m,z[[2,1]],z[[2,2]],z[[2,3]]]]<measured\[Mu]gLgR[[1,2]]&&measured\[Tau]gLgR[[2,1]]<gRSM+Re[\[Delta]gR[m,z[[3,1]],z[[3,2]],z[[3,3]]]]<measured\[Tau]gLgR[[2,2]]&&measured\[Mu]gLgR[[2,1]]<gRSM+Re[\[Delta]gR[m,z[[2,1]],z[[2,2]],z[[2,3]]]]<measured\[Mu]gLgR[[2,2]]

];


rD\[Mu]e=Compile[{{m,_Real},{cV,_Real},{cS,_Real},{cT,_Real}},
NIntegrate[bardhand\[CapitalGamma]BDl\[Nu]dqs[1+cV,runningS[m] cS,2runningT[m] cT,qs,2],{qs,m\[Mu]^2,(mB-mD)^2},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]/(NIntegrate[bardhand\[CapitalGamma]BDl\[Nu]dqs[1,0,0,qs,1],{qs,me^2,(mB-mD)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5]),
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True];

rD\[Mu]eLQ=Compile[{{m,_Real},{x,_Real,2},{y,_Real,2}},
Block[
{z,
cV\[Mu],cS\[Mu],cT\[Mu],
\[CapitalGamma]e,\[CapitalGamma]\[Mu]},

z=x.ConjugateTranspose[CKM];

cV\[Mu]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[2,2]]{x[[1,3]],x[[2,3]],x[[3,3]]};
cS\[Mu]=1/(2Sqrt[2]GFVcb) y[[2,2]]/(2m^2) {x[[1,3]],x[[2,3]],x[[3,3]]};
cT\[Mu]=-1/4 cS\[Mu];

\[CapitalGamma]\[Mu]=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDl\[Nu]dqs[KroneckerDelta[2,\[Nu]l]+cV\[Mu][[\[Nu]l]],runningS[m] cS\[Mu][[\[Nu]l]],2 runningT[m] cT\[Mu][[\[Nu]l]],qs,2],
{\[Nu]l,1,3}],
{qs,m\[Mu]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]e=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDl\[Nu]dqs[KroneckerDelta[1,\[Nu]l],0,0,qs,1],
{\[Nu]l,1,3}],
{qs,me^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]\[Mu]/\[CapitalGamma]e],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];

rD\[Mu]eLQCondition[m_,x_,y_]:=0.995-2Sqrt[0.022^2+0.039^2]<rD\[Mu]eLQ[m,x,y]<0.995+2Sqrt[0.022^2+0.039^2];


(* hep-ph/9806487 *)

R\[CapitalUpsilon]\[Nu]eCondition[m_,x_,y_]:=Block[{m\[CapitalUpsilon]=9.46,\[Alpha]=1/137,sw2=0.2223,g},
g=Sqrt[4\[Pi] \[Alpha]]/Sqrt[sw2];
(*\[CapitalGamma]\[CapitalUpsilon]ee=1.19 10^(-6);*)
((9GF^2 m\[CapitalUpsilon]^4)/(64\[Pi]^2 \[Alpha]^2) (2(-1+4/3 sw2)^2+(-1+4/3 sw2+2(mZ^2/m^2) (1-sw2)/g^2 Sum[Abs[x[[i,3]]]^2,{i,1,3}])^2))(*<((3 10^(-4))/((1.19 10^(-6))/(54.02 10^(-6))))*)];


(* From "Prospects of discovering new physics in rare charm decays" - Fajfer, Kosnik *)

brD\[Pi]\[Mu]\[Mu]Condition[m_,x_,y_]:=Block[{cll,crr,gll,grr,hll,hrr,
cc,aa,Br,NN,\[Tau]DGeV,\[Lambda]\[Lambda],T5,CT5,T,CPDash,CP,P,CSDash,CS,f0,S,C10,C10Dash,A,CT,C7,C7Dash,C9,C9Dash,
fPlus,fT,m\[Pi]=139.57 10^-3,fD,\[Tau]D=1040 10^-15,\[CapitalGamma]D,mc=0.54,\[Beta],z,V,\[Alpha]=1/127.},

fD=212 10^-3;
\[CapitalGamma]D=hbar/\[Tau]D;
\[Tau]DGeV=\[Tau]D/(6.58 10^-25);

z=x.ConjugateTranspose[CKM];

\[Beta][qq_]:=Sqrt[1-4m\[Mu]^2/qq];

fPlus[qq_]:=0.67/((1-qq/1.9^2)(1-0.28 qq/1.9^2));
fT[qq_]:=0.46/((1-qq/mD^2)(1-0.18 qq/mD^2));
f0[qq_]:=0.67/(1-1/1.27qq/1.9^2);

cll=-vev^2/(4m^2)z[[2,2]]Conjugate[z[[2,1]]];
crr=-vev^2/(4m^2)y[[2,2]]Conjugate[y[[2,1]]];
gll=runningSD[m] vev^2/(4m^2)Conjugate[y[[2,2]]]Conjugate[z[[2,1]]];
grr=runningSD[m] vev^2/(4m^2)z[[2,2]]y[[2,1]];
hll=-runningTD[m] gll/4;
hrr=-runningTD[m] grr/4;

C7=0;
C7Dash=0;
C9=(2\[Pi])/\[Alpha](cll);
C9Dash=(2\[Pi])/\[Alpha](crr);
CT=(2\[Pi])/\[Alpha](hll+hrr);
CT5=(2\[Pi])/\[Alpha](hll-hrr);
C10=(2\[Pi])/\[Alpha](-cll);
C10Dash=(2\[Pi])/\[Alpha](crr);
CS=(2\[Pi])/\[Alpha] gll;
CSDash=(2\[Pi])/\[Alpha] grr;
CP=(2\[Pi])/\[Alpha] gll;
CPDash=(-2\[Pi])/\[Alpha] grr;

\[Lambda]\[Lambda][qq_]:=(mD^2+m\[Pi]^2+qq)^2-4(mD^2m\[Pi]^2+m\[Pi]^2qq+mD^2qq);

V[qq_]:=(2 mc fT[qq])/(mD+m\[Pi])(C7+C7Dash)+fPlus[qq](C9 + C9Dash)+(8fT[qq]m\[Mu])/(mD+m\[Pi])CT;

A[qq_]:=fPlus[qq](C10+C10Dash);

S[qq_]:=(mD^2-m\[Pi]^2)/(2mc)f0[qq](CS+CSDash);

P[qq_]:=(mD^2-m\[Pi]^2)/(2mc)f0[qq](CP+CPDash)-m\[Mu](fPlus[qq]-(mD^2-m\[Pi]^2)/qq(f0[qq]-fPlus[qq]))(C10+C10Dash);

T[qq_]:=(2 fT[qq]\[Beta][qq] Sqrt[\[Lambda]\[Lambda][qq]] )/(mD+m\[Pi])CT;

T5[qq_]:=(2 fT[qq]\[Beta][qq] Sqrt[\[Lambda]\[Lambda][qq]] )/(mD+m\[Pi])CT5;

NN=(GF^2\[Alpha]^2)/((4\[Pi])^5mD^3); (*Took out \[Lambda]b because assume it will cancel with one from H*)

aa[qq_]:= \[Lambda]\[Lambda][qq]/2(Abs[V[qq]]^2+Abs[A[qq]]^2)+8m\[Mu]^2mD^2Abs[A[qq]]^2+2qq(\[Beta][qq]^2Abs[S[qq]]^2+Abs[P[qq]]^2)+4m\[Mu](mD^2-m\[Pi]^2+qq)Re[A[qq] Conjugate[P[qq]]] ;

cc[qq_]:=-((\[Lambda]\[Lambda][qq]\[Beta][qq]^2)/2)(Abs[V[qq]]^2+Abs[A[qq]]^2)+2qq(\[Beta][qq]^2Abs[T[qq]]^2+Abs[T5[qq]]^2)+4m\[Mu] \[Beta][qq] Sqrt[\[Lambda]\[Lambda][qq]]Re[V[qq]Conjugate[T[qq]]];

Br[qq_]:=\[Tau]DGeV 2 NN Sqrt[\[Lambda]\[Lambda][qq]]\[Beta][qq] (aa[qq]+cc[qq]/3);

(NIntegrate[Br[qq],{qq,1.56,2.99},WorkingPrecision->7,MinRecursion->4,AccuracyGoal->5]//Re )< 2.9 * 10^-8]



(* Grinstein *)

BrBc\[Tau]\[Nu][\[Epsilon]L_,\[Epsilon]P_]:=Module[{fB=434 10^-3,mB=6.275,md={0.005,0.095,4.18},mu={0.0022,0.54,160}},
(hbar/(0.507 10^-12))^-1 ((GF^2 mB Abs[CKM[[2,3]]]^2)/(8\[Pi]) fB^2 mLeps[[3]]^2 (1-mLeps[[3]]^2/mB^2)^2 Abs[1+\[Epsilon]L+\[Epsilon]P(mB^2/(mLeps[[3]](mu[[2]]+md[[3]])))]^2)
];

BrBc\[Tau]\[Nu]Running[\[Epsilon]L_,\[Epsilon]P_]:=BrBc\[Tau]\[Nu][\[Epsilon]L, 1.64769025413 \[Epsilon]P];

brBc\[Tau]\[Nu]Condition[m_,x_,y_]:=Block[{z,cV\[Tau],cS\[Tau]},
z=x.ConjugateTranspose[CKM];
cV\[Tau]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[3,2]]x[[;;,3]];
cS\[Tau]=1/(2Sqrt[2]GFVcb) y[[3,2]]/(2m^2) x[[;;,3]];
Sum[BrBc\[Tau]\[Nu][cV\[Tau][[i]],-cS\[Tau][[i]]],{i,3}]<0.3
];


(* Limit taken from recent Belle measurement Belleconf-1612 *)
rDStare\[Mu]LQCondition = Function[{m,x,y},
Block[
{z,
cV\[Mu],cS\[Mu],cT\[Mu],
\[CapitalGamma]e,\[CapitalGamma]\[Mu]},

z=x.ConjugateTranspose[CKM];

cV\[Mu]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[2,2]]x[[;;,3]];
cS\[Mu]=1/(2Sqrt[2]GFVcb) y[[2,2]]/(2m^2) x[[;;,3]];
cT\[Mu]=-1/4 cS\[Mu];

\[CapitalGamma]e=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDstarl\[Nu]dqs[KroneckerDelta[1,\[Nu]l],-(KroneckerDelta[1,\[Nu]l]),0,0,qs,1],
{\[Nu]l,1,3}],
{qs,mLeps[[1]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]\[Mu]=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDstarl\[Nu]dqs[KroneckerDelta[2,\[Nu]l]+cV\[Mu][[\[Nu]l]],-(cV\[Mu][[\[Nu]l]]+KroneckerDelta[2,\[Nu]l]),-runningS[m] cS\[Mu][[\[Nu]l]],2 runningT[m] cT\[Mu][[\[Nu]l]],qs,2],
{\[Nu]l,1,3}],
{qs,mLeps[[2]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5];

1.04-2Sqrt[0.05^2+0.01^2]<(\[CapitalGamma]e/\[CapitalGamma]\[Mu])<1.04+2Sqrt[0.05^2+0.01^2]]];


(* Implement NLO calculation from Michael's reference *)

(*brK\[Mu]\[Nu]Condition[m_,x_,y_]*)

(* eq 16 in 1605.07114 *)
rK\[Pi]l\[Nu][\[Epsilon]LD\[Mu]_,\[Epsilon]PD\[Mu]_]:=Block[{vusTilde,fK,f\[Pi],b0=2.5,fKonf\[Pi]Sq=1.194},
vusTilde=(1+\[Epsilon]LD\[Mu])CKM[[1,2]];
(Abs[vusTilde]^2 fKonf\[Pi]Sq)/Abs[CKM[[1,1]]]^2 (1-(2 b0)/m\[Mu] \[Epsilon]PD\[Mu])];

rK\[Pi]l\[Nu]OPS[cV_,cS_]:=rK\[Pi]l\[Nu][cV,-cS];

(* Need to find accurate experimental value! *)
rK\[Pi]l\[Nu]LQCondition[m_,x_,y_]:=Block[{cV,cS},
cV=(Conjugate[(x.ConjugateTranspose[CKM])[[2,1]]]x[[2,2]])/(2Sqrt[2]GF Abs[CKM[[1,2]]]2m^2);
cS=(y[[2,1]]x[[2,2]])/(2Sqrt[2]GF Abs[CKM[[1,2]]]2m^2);

0.6608 - 2 0.0029<rK\[Pi]l\[Nu]OPS[cV,cS]<0.6608 + 2 0.0029
];


myD00[m1s_,m2s_,m3s_,m4s_]:=With[{t1=m1s/m4s,t2=m2s/m4s,t3=m3s/m4s},
(-((t1^2 Log[t1])/(4(t1-1)(t1-t2)(t1-t3)))+(t2^2 Log[t2])/(4(t1-t2)(t2-1)(t2-t3))+(t3^2 Log[t3])/(4(t1-t3)(t3-1)(t3-t2))) 1/m4s
];
myD00[m1s_,m2s_,m3s_,m3s_]:=With[{t1=m1s/m3s,t2=m2s/m3s},
(-((t1^2 Log[t1])/((-1+t1)^2 (t1-t2)))+((1-t2)/(-1+t1)+(t2^2 Log[t2])/(t1-t2))/(-1+t2)^2)/(4 m3s)
];
myD00[m1s_,0,m3s_,m4s_]:=With[{t1=m1s/m3s,t3=m3s/m4s},
((t1-t1 t3) Log[t1]+(-1+t1) t3 Log[t3])/(4 m4s (-1+t1) (t1-t3) (-1+t3))
];
myD00[0,m3s_,m3s_,m4s_]:=With[{t3=m3s/m4s},
(1-t3+Log[t3])/(4 m4s (-1+t3)^2)
];
myD00[m2s_,m2s_,m3s_,m4s_]:=With[{t2=m2s/m4s,t3=m3s/m4s},
(t2 (-1+t3) (t2-2 t3+t2 t3) Log[t2]-(-1+t2) (t2 (t2-t3) (-1+t3)+(-1+t2) t3^2 Log[t3]))/(4 m4s (-1+t2)^2 (t2-t3)^2 (-1+t3))
];
myD00[m2s_,m2s_,m4s_,m4s_]:=With[{t2=m2s/m4s},
(1-t2^2+2 t2 Log[t2])/(4 m4s (-1+t2)^3)
];

Br\[Tau]\[Mu]\[Mu]\[Mu]=Function[{m\[Phi],x,y},Block[{t,\[Sigma]L,\[Sigma]R,A1L,A1R,A2L,A2R,B1L,B1R,B2L,B2R,B3L,B3R,FL,FR,\[CapitalGamma],FLL,FLR,FRR,FRL,ZL,ZR,a,\[Lambda]LeuMatrix,\[Lambda]ReuMatrix,e=Sqrt[4\[Pi] (1/137.)]},
\[Lambda]LeuMatrix=x.ConjugateTranspose[CKM];
\[Lambda]ReuMatrix=y;
Subscript[t, 1]=mu[[1]]^2/m\[Phi]^2;
Subscript[t, 2]=mu[[2]]^2/m\[Phi]^2;
Subscript[t, 3]=mu[[3]]^2/m\[Phi]^2;
A1L=Sum[(\[Lambda]LeuMatrix[[3,a]]Conjugate[\[Lambda]LeuMatrix[[2,a]]])/(384\[Pi]^2 m\[Phi]^2) 1/(Subscript[t, a]-1)^4 ((Subscript[t, a]-1)(10+(-17+Subscript[t, a])Subscript[t, a])-2(4-6Subscript[t, a]-Subscript[t, a]^3)Log[Subscript[t, a]]),{a,1,3}];
A1R=Sum[(\[Lambda]ReuMatrix[[3,a]]Conjugate[\[Lambda]ReuMatrix[[2,a]]])/(384\[Pi]^2 m\[Phi]^2) 1/(Subscript[t, a]-1)^4 ((Subscript[t, a]-1)(10+(-17+Subscript[t, a])Subscript[t, a])-2(4-6Subscript[t, a]-Subscript[t, a]^3)Log[Subscript[t, a]]),{a,1,3}];
\[Sigma]L=Sum[1/(16\[Pi]^2) 1/m\[Phi]^2 \[Lambda]LeuMatrix[[3,a]]Conjugate[\[Lambda]LeuMatrix[[2,a]]]m\[Tau] (1+4Subscript[t, a]-5Subscript[t, a]^2+2Subscript[t, a](2+Subscript[t, a])Log[Subscript[t, a]])/(4(Subscript[t, a]-1)^4)+1/(16\[Pi]^2) 1/m\[Phi]^2 \[Lambda]LeuMatrix[[3,a]]Conjugate[\[Lambda]ReuMatrix[[2,a]]]mu[[a]] (7-8Subscript[t, a]+Subscript[t, a]^2+2(2+Subscript[t, a])Log[Subscript[t, a]])/(2(Subscript[t, a]-1)^3),{a,1,3}];
\[Sigma]R=Sum[1/(16\[Pi]^2) 1/m\[Phi]^2 \[Lambda]ReuMatrix[[3,a]]Conjugate[\[Lambda]ReuMatrix[[2,a]]]m\[Tau] (1+4Subscript[t, a]-5Subscript[t, a]^2+2Subscript[t, a](2+Subscript[t, a])Log[Subscript[t, a]])/(4(Subscript[t, a]-1)^4)+1/(16\[Pi]^2) 1/m\[Phi]^2 \[Lambda]ReuMatrix[[3,a]]Conjugate[\[Lambda]LeuMatrix[[2,a]]]mu[[a]] (7-8Subscript[t, a]+Subscript[t, a]^2+2(2+Subscript[t, a])Log[Subscript[t, a]])/(2(Subscript[t, a]-1)^3),{a,1,3}];
A2R=\[Sigma]L/m\[Tau];
A2L=\[Sigma]R/m\[Tau];
FL=Sum[- ((3e \[Lambda]LeuMatrix[[3,a]]Conjugate[\[Lambda]LeuMatrix[[2,a]]])/(32\[Pi]^3 Sin[\[Theta]w]Cos[\[Theta]w])) (Subscript[t, a](1-Subscript[t, a]+Log[Subscript[t, a]]))/(Subscript[t, a]-1)^2,{a,1,3}];
FR=Sum[- ((3e \[Lambda]ReuMatrix[[3,a]]Conjugate[\[Lambda]ReuMatrix[[2,a]]])/(32\[Pi]^3 Sin[\[Theta]w]Cos[\[Theta]w])) (Subscript[t, a](1-Subscript[t, a]+Log[Subscript[t, a]]))/(Subscript[t, a]-1)^2,{a,1,3}];
ZL=-e/(Sin[\[Theta]w]Cos[\[Theta]w]) (-1/2+Sin[\[Theta]w]^2);
ZR=-e/(Sin[\[Theta]w]Cos[\[Theta]w]) Sin[\[Theta]w]^2;
FLL=(FL ZL)/(e^2 mZ^2);
FRR=(FR ZR)/(e^2 mZ^2);
FLR=(FL ZR)/(e^2 mZ^2);
FRL=(FR ZL)/(e^2 mZ^2);
B1L=Sum[-(3/(16\[Pi]^2 e^2))\[Lambda]LeuMatrix[[3,i]]Conjugate[\[Lambda]LeuMatrix[[2,i]]]Conjugate[\[Lambda]LeuMatrix[[2,j]]]\[Lambda]LeuMatrix[[2,j]] myD00[m\[Phi]^2,m\[Phi]^2,mu[[i]]^2,mu[[j]]^2],{i,1,3},{j,1,3}];
B1R=Sum[-(3/(16\[Pi]^2 e^2))\[Lambda]ReuMatrix[[3,i]]Conjugate[\[Lambda]ReuMatrix[[2,i]]]Conjugate[\[Lambda]ReuMatrix[[2,j]]]\[Lambda]ReuMatrix[[2,j]] myD00[m\[Phi]^2,m\[Phi]^2,mu[[i]]^2,mu[[j]]^2],{i,1,3},{j,1,3}];
B2L=Sum[-(3/(16\[Pi]^2 e^2))\[Lambda]LeuMatrix[[3,i]]Conjugate[\[Lambda]LeuMatrix[[2,i]]]Conjugate[\[Lambda]ReuMatrix[[2,j]]]\[Lambda]ReuMatrix[[2,j]] myD00[m\[Phi]^2,m\[Phi]^2,mu[[i]]^2,mu[[j]]^2],{i,1,3},{j,1,3}];
B2R=Sum[-(3/(16\[Pi]^2 e^2))\[Lambda]ReuMatrix[[3,i]]Conjugate[\[Lambda]ReuMatrix[[2,i]]]Conjugate[\[Lambda]LeuMatrix[[2,j]]]\[Lambda]LeuMatrix[[2,j]] myD00[m\[Phi]^2,m\[Phi]^2,mu[[i]]^2,mu[[j]]^2],{i,1,3},{j,1,3}];

\[CapitalGamma]=hbar/(290 10^-15);
1/\[CapitalGamma] e^4/(512\[Pi]^3) m\[Tau]^5 (Abs[A1L]^2+Abs[A1R]^2-2(A1L Conjugate[A2R]+A2L Conjugate[A1R] + Conjugate[A1L Conjugate[A2R]+A2L Conjugate[A1R]]) +(Abs[A2L]^2+Abs[A2R]^2)(16/3 Log[m\[Tau]/m\[Mu]]-22/3)+1/6 (Abs[B1L]^2+Abs[B1R]^2)+1/3 (Abs[B2L]^2+Abs[B2R]^2)+1/3 (A1L Conjugate[B1L]+A1R Conjugate[B1R]+A1L Conjugate[B2L]+A1R Conjugate[B2R] + Conjugate[A1L Conjugate[B1L]+A1R Conjugate[B1R]+A1L Conjugate[B2L]+A1R Conjugate[B2R]])-2/3 (A2R Conjugate[B1L]+A2L Conjugate[B1R]+A2L Conjugate[B2R]+A2R Conjugate[B2L]+Conjugate[A2R Conjugate[B1L]+A2L Conjugate[B1R]+A2L Conjugate[B2R]+A2R Conjugate[B2L]])+1/3 (2(Abs[FLL]^2+Abs[FRR]^2)+Abs[FLR]^2+Abs[FRL]^2+(B1L Conjugate[FLL]+B1R Conjugate[FRR] + B2L Conjugate[FLR]+B2R Conjugate[FRL]+Conjugate[B1L Conjugate[FLL]+B1R Conjugate[FRR] + B2L Conjugate[FLR]+B2R Conjugate[FRL]])+2(A1L Conjugate[FLL]+A1R Conjugate[FRR]+Conjugate[A1L Conjugate[FLL]+A1R Conjugate[FRR]])+(A1L Conjugate[FLR]+A1R Conjugate[FRL]+Conjugate[A1L Conjugate[FLR]+A1R Conjugate[FRL]]) - 4(A2R Conjugate[FLL]+A2L Conjugate[FRR]+Conjugate[A2R Conjugate[FLL]+A2L Conjugate[FRR]])-2(A2L Conjugate[FRL]+A2L Conjugate[FLR]+Conjugate[A2L Conjugate[FRL]+A2L Conjugate[FLR]])))]];

Br\[Tau]\[Mu]\[Mu]\[Mu]Condition[m_,x_,y_]:=Re[Br\[Tau]\[Mu]\[Mu]\[Mu][m,x,y]]<2.1 10^-8;

BN\[Tau]\[Mu]\[Gamma]Condition[m\[Phi]_,x_,\[Lambda]ReuMatrix_]:=With[{ac=1+0.17Log[m\[Phi]/1000],at=1+1.06Log[m\[Phi]/1000],\[Lambda]LeuMatrix=x.ConjugateTranspose[CKM]},
\[Sqrt]Abs[ac \[Lambda]ReuMatrix[[3,2]]Conjugate[\[Lambda]LeuMatrix[[2,2]]]+ 20.7 at \[Lambda]ReuMatrix[[3,3]]Conjugate[\[Lambda]LeuMatrix[[2,3]]] - 0.015 (Sum[\[Lambda]LeuMatrix[[2,i]]Conjugate[\[Lambda]LeuMatrix[[3,i]]],{i,3}])]^2+Abs[ac \[Lambda]LeuMatrix[[3,2]]Conjugate[\[Lambda]ReuMatrix[[2,2]]]+ 20.7 at \[Lambda]LeuMatrix[[3,3]]Conjugate[\[Lambda]ReuMatrix[[2,3]]] - 0.015 (Sum[\[Lambda]ReuMatrix[[2,i]]Conjugate[\[Lambda]ReuMatrix[[3,i]]],{i,3}])]^2<0.017 (m\[Phi]/1000)^2
];


vanillaPassConstraintsQ=Function[{m,x,y},
rK\[Nu]\[Nu]Condition[m,x]&&
cBsCondition[m,x]&&
brD\[Mu]\[Mu]Condition[m,x,y]&&
brB\[Tau]\[Nu]Condition[m,x,y]&&
brB\[Mu]\[Nu]Condition[m,x,y]&&
brK\[Mu]\[Nu]Condition[m,x,y]&&
rKe\[Mu]Condition[m,x,y]&&
rK\[Tau]\[Mu]Condition[m,x,y]&&
brDs\[Mu]\[Nu]Condition[m,x,y]&&
brDs\[Tau]\[Nu]Condition[m,x,y]&&
zllConditionsEWFit[m,x,y]&&
rD\[Mu]eLQCondition[m,x,y]&&
brD\[Pi]\[Mu]\[Mu]Condition[m,x,y]&&
brBc\[Tau]\[Nu]Condition[m,x,y]&&
rDStare\[Mu]LQCondition[m,x,y]&&
Br\[Tau]\[Mu]\[Mu]\[Mu]Condition[m,x,y]&&
BN\[Tau]\[Mu]\[Gamma]Condition[m,x,y]
];
