(* ::Package:: *)

Print["BAnomalies Version: 1.0 (8 Feb 2017)."];
Print["Author: John Gargalionis"];
Print["Please cite the papers whose results we use:\n -- Bardhan, Byakti, Ghosh, 'A closer look at the RD and RD* anomalies (arXiv:1610.03038v2)'\n -- Becirevic, Kosnik, Sumensari, Funchal, 'Palatable leptoquark scenarios for lepton flavor violation in Exclusive b -> s l_1 l_2 decays (arXiv:1608.07583v2)'\n and references therein."];
Print["Note that the package BAnomalies needs to be loaded after PhysicsWorkbench."];


BeginPackage["BAnomalies`"];

{bardhand\[CapitalGamma]BDl\[Nu]dqs, rD, rDLQ, bardhand\[CapitalGamma]BDstarl\[Nu]dqs, rDStar, rDStarLQ, CLL, CLR, dBrdqsRK, rKLQ}

EndPackage[];


<<"~/Dropbox/mathematicaPackages/RunDec.m"

\[Alpha]sMt5 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mt/.NumDef,5,3];
\[Alpha]sMt6 = AlphasExact[asMz/.NumDef,Mz/.NumDef,Mt/.NumDef,6,3];
\[Alpha]sMb = AlphasExact[asMz/.NumDef,Mz/.NumDef,4.2,5,3];
\[Alpha]sMLQ[m_] := AsRunDec[asMz/.NumDef,Mz/.NumDef,m,3];

(* running factors for scalar and tensor operators *)
runningS[m_]:=(\[Alpha]sMt5/\[Alpha]sMb)^(-12/23)(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(-12/21);
runningT[m_]:=(\[Alpha]sMt5/\[Alpha]sMb)^(4/23)(\[Alpha]sMLQ[m]/\[Alpha]sMt6)^(4/21);

lambda[a_,b_,c_]:=a^2+b^2+c^2-2(a b + b c + a c);


bardhand\[CapitalGamma]BDl\[Nu]dqsMinus=Function[{cvl,csl,ctl,qs,l},
Block[
{magpD,
alDMinus,clDMinus,alDPlus,clDPlus,
z,r,\[Phi]Plus,\[Phi]0,aPlus,a0,
fPlus,fT,f0,
nPrefactor,branchingRatioMinus,branchingRatioPlus},
magpD=Sqrt[lambda[mB^2,mD^2,qs]]/(2mB);

z=(Sqrt[(mB+mD)^2-qs]-2Sqrt[mB mD])/(Sqrt[(mB+mD)^2-qs]+2Sqrt[mB mD]);
r=mD/mB;
\[Phi]Plus=1.1213 (1+z)^2/((1+r)(1-z)+2Sqrt[r](1+z))^5;
\[Phi]0=0.5299 ((1+z)(1-z)^(3/2))/((1+r)(1-z)+2Sqrt[r](1+z))^4;
Subscript[aPlus, 0]=0.01261;Subscript[aPlus, 1]=-0.0963;Subscript[aPlus, 2]=0.37;Subscript[aPlus, 3]=-0.05;

Subscript[a0, 0]=0.01140;
Subscript[a0, 1]=-0.0590;
Subscript[a0, 2]=0.19;
Subscript[a0, 3]=-0.03;

fPlus=1/\[Phi]Plus Sum[Subscript[aPlus, k] z^k,{k,0,3}];
f0=1/\[Phi]0 Sum[Subscript[a0, k] z^k,{k,0,3}];
fT=0.69/((1-qs/6.4^2)(1-0.56 qs/6.4^2)); 

alDMinus=8((mB^2 magpD^2)/qs Abs[cvl]^2 fPlus^2+mLeps[[l]] (4mB^2 magpD^2)/(qs (mB + mD)) Re[ctl Conjugate[cvl] fPlus fT]+mLeps[[l]]^2 (4magpD^2 mB^2)/(qs (mB+mD)^2) Abs[ctl]^2 fT^2);
clDMinus=8(-((mB^2 magpD^2)/qs) Abs[cvl]^2 fPlus^2-mLeps[[l]] (4magpD^2 mB^2)/(qs(mB+mD)) Re[cvl ctl]fPlus fT-mLeps[[l]]^2 (4magpD^2 mB^2)/((mB+mD)^2 qs) Abs[ctl]^2 fT^2);

\[Tau]B=(1.638 10^-12)/(6.58 10^-25); (* lifetime of B meson in GeV^-1 *)

nPrefactor=(\[Tau]B Abs[GFVcb]^2 qs)/(256 \[Pi]^3 mB^2) (1-mLeps[[l]]^2/qs)^2;
branchingRatioMinus=nPrefactor magpD (2 alDMinus + 2/3 clDMinus);

branchingRatioMinus

]
];

bardhand\[CapitalGamma]BDl\[Nu]dqsPlus=Function[{cvl,csl,ctl,qs,l},
Block[
{magpD,
alDMinus,clDMinus,alDPlus,clDPlus,
z,r,\[Phi]Plus,\[Phi]0,aPlus,a0,
fPlus,fT,f0,
nPrefactor,branchingRatioMinus,branchingRatioPlus},
magpD=Sqrt[lambda[mB^2,mD^2,qs]]/(2mB);

z=(Sqrt[(mB+mD)^2-qs]-2Sqrt[mB mD])/(Sqrt[(mB+mD)^2-qs]+2Sqrt[mB mD]);
r=mD/mB;
\[Phi]Plus=1.1213 (1+z)^2/((1+r)(1-z)+2Sqrt[r](1+z))^5;
\[Phi]0=0.5299 ((1+z)(1-z)^(3/2))/((1+r)(1-z)+2Sqrt[r](1+z))^4;
Subscript[aPlus, 0]=0.01261;Subscript[aPlus, 1]=-0.0963;Subscript[aPlus, 2]=0.37;Subscript[aPlus, 3]=-0.05;

Subscript[a0, 0]=0.01140;
Subscript[a0, 1]=-0.0590;
Subscript[a0, 2]=0.19;
Subscript[a0, 3]=-0.03;

fPlus=1/\[Phi]Plus Sum[Subscript[aPlus, k] z^k,{k,0,3}];
f0=1/\[Phi]0 Sum[Subscript[a0, k] z^k,{k,0,3}];
fT=0.69/((1-qs/6.4^2)(1-0.56 qs/6.4^2)); 

alDPlus=8((mB^2-mD^2)^2/(4(mb-mc)^2) Abs[csl]^2 f0^2+mLeps[[l]] (mB^2-mD^2)^2/(2qs (mb-mc)) Re[csl Conjugate[cvl]]f0^2+mLeps[[l]]^2 (mB^2-mD^2)^2/(4qs^2) Abs[cvl]^2 f0^2);

clDPlus=8((4 magpD^2 mB^2)/(mB+mD)^2 Abs[ctl]^2 fT^2+mLeps[[l]] (4mB^2 magpD^2)/((mB+mD)qs) Re[cvl Conjugate[ctl]]fPlus fT+mLeps[[l]]^2 (magpD^2 mB^2)/qs^2 Abs[cvl]^2 fPlus^2);

\[Tau]B=(1.638 10^-12)/(6.58 10^-25); (* lifetime of B meson in GeV^-1 *)

nPrefactor=(\[Tau]B Abs[GFVcb]^2 qs)/(256 \[Pi]^3 mB^2) (1-mLeps[[l]]^2/qs)^2;
branchingRatioPlus=nPrefactor magpD (2 alDPlus + 2/3 clDPlus);

branchingRatioPlus

]
];


bardhand\[CapitalGamma]BDl\[Nu]dqs=Function[{cvl,csl,ctl,qs,l},
Block[
{magpD,
alDMinus,clDMinus,alDPlus,clDPlus,
z,r,\[Phi]Plus,\[Phi]0,aPlus,a0,
fPlus,fT,f0,
nPrefactor,branchingRatioMinus,branchingRatioPlus},
magpD=Sqrt[lambda[mB^2,mD^2,qs]]/(2mB);

z=(Sqrt[(mB+mD)^2-qs]-2Sqrt[mB mD])/(Sqrt[(mB+mD)^2-qs]+2Sqrt[mB mD]);
r=mD/mB;
\[Phi]Plus=1.1213 (1+z)^2/((1+r)(1-z)+2Sqrt[r](1+z))^5;
\[Phi]0=0.5299 ((1+z)(1-z)^(3/2))/((1+r)(1-z)+2Sqrt[r](1+z))^4;
Subscript[aPlus, 0]=0.01261;Subscript[aPlus, 1]=-0.0963;Subscript[aPlus, 2]=0.37;Subscript[aPlus, 3]=-0.05;

Subscript[a0, 0]=0.01140;
Subscript[a0, 1]=-0.0590;
Subscript[a0, 2]=0.19;
Subscript[a0, 3]=-0.03;

fPlus=1/\[Phi]Plus Sum[Subscript[aPlus, k] z^k,{k,0,3}];
f0=1/\[Phi]0 Sum[Subscript[a0, k] z^k,{k,0,3}];
fT=0.69/((1-qs/6.4^2)(1-0.56 qs/6.4^2)); 

alDMinus=8((mB^2 magpD^2)/qs Abs[cvl]^2 fPlus^2+mLeps[[l]] (4mB^2 magpD^2)/(qs (mB + mD)) Re[ctl Conjugate[cvl] fPlus fT]+mLeps[[l]]^2 (4magpD^2 mB^2)/(qs (mB+mD)^2) Abs[ctl]^2 fT^2);
clDMinus=8(-((mB^2 magpD^2)/qs) Abs[cvl]^2 fPlus^2-mLeps[[l]] (4magpD^2 mB^2)/(qs(mB+mD)) Re[cvl ctl]fPlus fT-mLeps[[l]]^2 (4magpD^2 mB^2)/((mB+mD)^2 qs) Abs[ctl]^2 fT^2);

alDPlus=8((mB^2-mD^2)^2/(4(mb-mc)^2) Abs[csl]^2 f0^2+mLeps[[l]] (mB^2-mD^2)^2/(2qs (mb-mc)) Re[csl Conjugate[cvl]]f0^2+mLeps[[l]]^2 (mB^2-mD^2)^2/(4qs^2) Abs[cvl]^2 f0^2);

clDPlus=8((4 magpD^2 mB^2)/(mB+mD)^2 Abs[ctl]^2 fT^2+mLeps[[l]] (4mB^2 magpD^2)/((mB+mD)qs) Re[cvl Conjugate[ctl]]fPlus fT+mLeps[[l]]^2 (magpD^2 mB^2)/qs^2 Abs[cvl]^2 fPlus^2);

\[Tau]B=(1.638 10^-12)/(6.58 10^-25); (* lifetime of B meson in GeV^-1 *)

nPrefactor=(\[Tau]B Abs[GFVcb]^2 qs)/(256 \[Pi]^3 mB^2) (1-mLeps[[l]]^2/qs)^2;
branchingRatioMinus=nPrefactor magpD (2 alDMinus + 2/3 clDMinus);
branchingRatioPlus=nPrefactor magpD (2 alDPlus + 2/3 clDPlus);

branchingRatioMinus+branchingRatioPlus

]
];

rDLQ=Compile[{{m,_Real},{x,_Real,2},{y,_Real,2}},
Block[
{z,
cV\[Tau],cS\[Tau],cT\[Tau],
cV\[Mu],cS\[Mu],cT\[Mu],
\[CapitalGamma]\[Tau],\[CapitalGamma]\[Mu]},

z=x.ConjugateTranspose[CKM];

cV\[Tau]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[3,2]]{x[[1,3]],x[[2,3]],x[[3,3]]};
cS\[Tau]=1/(2Sqrt[2]GFVcb) y[[3,2]]/(2m^2) {x[[1,3]],x[[2,3]],x[[3,3]]};
cT\[Tau]=(-1/4) cS\[Tau];

cV\[Mu]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[2,2]]{x[[1,3]],x[[2,3]],x[[3,3]]};
cS\[Mu]=1/(2Sqrt[2]GFVcb) y[[2,2]]/(2m^2) {x[[1,3]],x[[2,3]],x[[3,3]]};
cT\[Mu]=(-1/4) cS\[Mu];

\[CapitalGamma]\[Tau]=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDl\[Nu]dqs[KroneckerDelta[3,\[Nu]l]+cV\[Tau][[\[Nu]l]],runningS[m] cS\[Tau][[\[Nu]l]],2 runningT[m] cT\[Tau][[\[Nu]l]],qs,3],
{\[Nu]l,1,3}],
{qs,mLeps[[3]]^2,(mB-mDStar)^2},MinRecursion->9,WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]\[Mu]=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDl\[Nu]dqs[KroneckerDelta[2,\[Nu]l]+cV\[Mu][[\[Nu]l]],runningS[m] cS\[Mu][[\[Nu]l]],2 runningT[m] cT\[Mu][[\[Nu]l]],qs,2],
{\[Nu]l,1,3}],
{qs,mLeps[[2]]^2,(mB-mDStar)^2},MinRecursion->9,WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]\[Tau]/\[CapitalGamma]\[Mu]],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];

rD=Compile[{{m,_Real},{cV,_Real},{cS,_Real},{cT,_Real}},
Block[{},
NIntegrate[bardhand\[CapitalGamma]BDl\[Nu]dqs[1+cV,runningS[m] cS,2runningT[m] cT,qs,3],{qs,mLeps[[3]]^2,(mB-mD)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5]/(NIntegrate[bardhand\[CapitalGamma]BDl\[Nu]dqs[1,0,0,qs,2],{qs,mLeps[[2]]^2,(mB-mD)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5])],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True];

rD\[Mu]mLQ=Compile[{{cV,_Real},{cS,_Real},{cT,_Real}},
Block[{runningS,runningT},
runningS=1;
runningT=1;
NIntegrate[bardhand\[CapitalGamma]BDl\[Nu]dqs[1+cV,runningS cS,2runningT cT,qs,3],{qs,mLeps[[3]]^2,(mB-mD)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5]/(NIntegrate[bardhand\[CapitalGamma]BDl\[Nu]dqs[1,0,0,qs,2],{qs,mLeps[[2]]^2,(mB-mD)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]9,*)WorkingPrecision->7,AccuracyGoal->5])],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True];


bardhand\[CapitalGamma]BDstarl\[Nu]dqs=Function[{cvl,cal,cpl,ctl,qs,l},
Block[
{magpDstar,\[Lambda],
r,w,z,
r11,r21,\[Rho]DstarSq,ha11,
ha1,r1,r2,r3,
hv,ha2,ha3,ht1,ht2,ht3,
v,a0,a1,a2,t1,t2,t3,
alDstarMinus,clDstarMinus,alDstarPlus,clDstarPlus,
nPrefactor,branchingRatioMinus,branchingRatioPlus},

magpDstar=Sqrt[lambda[mB^2,mDStar^2,qs]]/(2mB);
(* Not sure if this is correct *)
\[Lambda]=lambda[mB^2,mDStar^2,qs];

r=mDStar/mB;
w=(mB^2+mDStar^2-qs)/(2mB mDStar);
z=(Sqrt[w+1]-Sqrt[2])/(Sqrt[w+1]+Sqrt[2]);

r11=1.406;
r21=0.853;
\[Rho]DstarSq=1.207;
ha11=0.906;

ha1=ha11 (1-8 \[Rho]DstarSq z + (53\[Rho]DstarSq - 15)z^2-(231\[Rho]DstarSq - 91)z^3);
r1=r11-0.12(w-1)+0.05(w-1)^2;
r2=r21+0.11(w-1)-0.06(w-1)^2;
r3=1.22 -0.052(w-1)+0.026(w-1)^2;

hv=r1 ha1;
ha2=(r2-r3)/(2r) ha1;
ha3=(r2+r3)/2 ha1;
ht1=1/(2(1+r^2-2r w)) ((mb-mc)/(mB-mDStar) (1-r)^2 (w+1)ha1-(mb+mc)/(mB+mDStar) (1+r)^2 (w-1)hv);
ht2=((1-r^2)(w+1))/(2(1+r^2-2r w)) ((mb-mc)/(mB-mDStar) ha1-(mb+mc)/(mB+mDStar) hv);
ht3=-(1/(2(1+r)(1+r^2-2r w)))(2 (mb-mc)/(mB-mDStar) r(w+1)ha1-(mb-mc)/(mB-mDStar) (1+r^2-2r w)(ha3-r ha2)-(mb+mc)/(mB+mDStar) (1+r)^2 hv);

v=(mB+mDStar)/(2Sqrt[mB mDStar]) hv;
a1=((mB+mDStar)^2-qs)/(2Sqrt[mB mDStar](mB+mDStar)) ha1;
a2=(mB+mDStar)/(2Sqrt[mB mDStar]) (ha3+mDStar/mB ha2);
a0=1/(2Sqrt[mB mDStar]) (((mB+mDStar)^2-qs)/(2mDStar) ha1-(mB^2-mDStar^2+qs)/(2mB) ha2-(mB^2-mDStar^2-qs)/(2mDStar) ha3);
t1=1/(2Sqrt[mB mDStar]) ((mB+mDStar)ht1-(mB-mDStar)ht2);
t2=1/(2Sqrt[mB mDStar]) (((mB+mDStar)^2-qs)/(mB+mDStar) ht1-((mB-mDStar)^2-qs)/(mB-mDStar) ht2);
t3=1/(2Sqrt[mB mDStar]) ((mB-mDStar)ht1-(mB+mDStar)ht2-(2(mB^2-mDStar^2))/mB ht3);

alDstarMinus=(8mB^2 magpDstar^2)/(mB+mDStar)^2 Abs[cvl]^2 v^2+((mB+mDStar)^2 (8mDStar^2 qs + \[Lambda]))/(2mDStar^2 qs) Abs[cal]^2 a1^2+(8mB^4 magpDstar^4)/(mDStar^2 (mB+mDStar)^2 qs) Abs[cal]^2 a2^2-(4magpDstar^2 mB^2 (mB^2-mDStar^2-qs))/(mDStar^2 qs) Abs[cal]^2 a1 a2+mLeps[[l]]((32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1+(8(mB+mDStar)(2mDStar^2 (mB^2-mDStar^2)+mB^2 magpDstar^2))/(qs mDStar^2) Re[cal Conjugate[ctl]]a1 t2-(8mB^2 (mB^2-mDStar^2-qs)magpDstar^2)/(qs (mB-mDStar)mDStar^2) Re[cal Conjugate[ctl]]a1 t3-(8mB^2 (mB^2+3mDStar^2-qs)magpDstar^2)/(qs(mB+mDStar)mDStar^2) Re[cal Conjugate[ctl]]a2 t2+(32mB^4 magpDstar^4)/(qs mDStar^2 (mB+mDStar)(mB^2-mDStar^2)) Re[cal Conjugate[ctl]]a2 t3)+mLeps[[l]]^2 ((32mB^2 magpDstar^2)/qs^2 Abs[ctl]^2 t1^2+(2(8mDStar^2 (2(mB^2+mDStar^2)-qs)qs+(4mDStar^2+qs)\[Lambda]))/(qs^2 mDStar^2) Abs[ctl]^2 t2^2+(32mB^4 magpDstar^4)/(qs mDStar^2 (mB^2-mDStar^2)^2) Abs[ctl]^2 t3^2-(16mB^2 magpDstar^2 (mB^2+3mDStar^2-qs))/(qs mDStar^2 (mB^2-mDStar^2)) Abs[ctl]^2 t2 t3);

clDstarMinus=(8magpDstar^2 mB^2)/(mB+mDStar)^2 Abs[cvl]^2 v^2-((mB+mDStar)^2 \[Lambda])/(2mDStar^2 qs) Abs[cal]^2 a1^2-(8 magpDstar^4 mB^4)/((mB+mDStar)^2 mDStar^2 qs) Abs[cal]^2 a2^2+(4magpDstar^2 mB^2 (mB^2-mDStar^2-qs))/(mDStar^2 qs) Abs[cal]^2 a1 a2+mLeps[[l]]((32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1-(8mB^2 (mB+mDStar)magpDstar^2)/(qs mDStar^2) Re[cal Conjugate[ctl]]a1 t2+(8mB^2 (mB^2-mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB-mDStar)) Re[cal Conjugate[ctl]]a1 t3 + (8mB^2 (mB^2+3mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB+mDStar)) Re[cal Conjugate[ctl]]a2 t2-(32mB^4 magpDstar^4)/(qs mDStar^2 (mB+mDStar)(mB^2-mDStar^2)) Re[cal Conjugate[ctl]]a2 t3)+mLeps[[l]]^2 ((32mB^2 magpDstar^2)/qs^2 Abs[ctl]^2 t1^2+(2(4mDStar^2-qs)\[Lambda])/(mDStar^2 qs^2) Abs[ctl]^2 t2^2-(32mB^4 magpDstar^4)/(qs mDStar^2 (mB^2-mDStar^2)^2) Abs[ctl]^2 t3^2+(16mB^2 magpDstar^2 (mB^2+3mDStar^2-qs))/(qs mDStar^2 (mB^2-mDStar^2)) Abs[ctl]^2 t2 t3);

alDstarPlus=(8magpDstar^2 mB^2)/(mb+mc)^2 Abs[cpl]^2 a0^2+(32mB^2 magpDstar^2)/qs Abs[ctl]^2 t1^2+(8(mB^2-mDStar^2)^2)/qs Abs[ctl]^2 t2^2-mLeps[[l]]((16 magpDstar^2 mB^2)/((mb+mc)qs) Re[cal Conjugate[cpl]]a0^2-(32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1-(8(mB+mDStar)(mB^2-mDStar^2))/qs Re[cal Conjugate[ctl]]a1 t2)+mLeps[[l]]^2 ((8magpDstar^2 mB^2)/qs^2 Abs[cal]^2 a0^2+(8magpDstar^2 mB^2)/((mB+mDStar)^2 qs) Abs[cvl]^2 v^2+(2(mB+mDStar)^2)/qs Abs[cal]^2 a1^2);

clDstarPlus=-((32mB^2 magpDstar^2)/qs) Abs[ctl]^2 t1^2-(2(4mDStar^2-qs)\[Lambda])/(mDStar^2 qs) Abs[ctl]^2 t2^2+(32mB^4 magpDstar^4)/(mDStar^2 (mB^2-mDStar^2)^2) Abs[ctl]^2 t3^2-(16mB^2 magpDstar^2 (mB^2+3mDStar^2-qs))/(mDStar^2 (mB^2-mDStar^2)) Abs[ctl]^2 t2 t3-mLeps[[l]]((32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1 - (8mB^2 (mB+mDStar)magpDstar^2)/(qs mDStar^2) Re[cal Conjugate[ctl]]a1 t2+(8mB^2 (mB^2-mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB-mDStar)) Re[cal Conjugate[ctl]]a1 t3+(8mB^2 (mB^2+3mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB+mDStar)) Re[cal Conjugate[ctl]]a2 t2-(32mB^4 magpDstar^4)/(qs mDStar^2 (mB+mDStar)(mB^2-mDStar^2)) Re[cal Conjugate[ctl]]a2 t3)+mLeps[[l]]^2 (-((8magpDstar^2 mB^2)/((mB+mDStar)^2 qs)) Abs[cvl]^2 v^2+((mB+mDStar)^2 \[Lambda])/(2mDStar^2 qs^2) Abs[cal]^2 a1^2+(8magpDstar^4 mB^4)/(mDStar^2 (mB+mDStar)^2 qs^2) Abs[cal]^2 a2^2-(4magpDstar^2 mB^2)/(mDStar^2 qs^2) (mB^2-mDStar^2-qs)Abs[cal]^2 a1 a2);

nPrefactor=(\[Tau]B Abs[GFVcb]^2 qs)/(256 \[Pi]^3 mB^2) (1-mLeps[[l]]^2/qs)^2;
branchingRatioMinus=nPrefactor magpDstar (2 alDstarMinus + 2/3 clDstarMinus);
branchingRatioPlus=nPrefactor magpDstar (2 alDstarPlus + 2/3 clDstarPlus);

branchingRatioMinus+branchingRatioPlus
]
];

bardhand\[CapitalGamma]BDstarl\[Nu]dqsMinus=Function[{cvl,cal,cpl,ctl,qs,l},
Block[
{magpDstar,\[Lambda],
r,w,z,
r11,r21,\[Rho]DstarSq,ha11,
ha1,r1,r2,r3,
hv,ha2,ha3,ht1,ht2,ht3,
v,a0,a1,a2,t1,t2,t3,
alDstarMinus,clDstarMinus,alDstarPlus,clDstarPlus,
nPrefactor,branchingRatioMinus,branchingRatioPlus,mb=4.18,mc=1.275},

magpDstar=Sqrt[lambda[mB^2,mDStar^2,qs]]/(2mB);
(* Not sure if this is correct *)
\[Lambda]=lambda[mB^2,mDStar^2,qs];

r=mDStar/mB;
w=(mB^2+mDStar^2-qs)/(2mB mDStar);
z=(Sqrt[w+1]-Sqrt[2])/(Sqrt[w+1]+Sqrt[2]);

r11=1.406;
r21=0.853;
\[Rho]DstarSq=1.207;
ha11=0.906;

ha1=ha11 (1-8 \[Rho]DstarSq z + (53\[Rho]DstarSq - 15)z^2-(231\[Rho]DstarSq - 91)z^3);
r1=r11-0.12(w-1)+0.05(w-1)^2;
r2=r21+0.11(w-1)-0.06(w-1)^2;
r3=1.22 -0.052(w-1)+0.026(w-1)^2;

hv=r1 ha1;
ha2=(r2-r3)/(2r) ha1;
ha3=(r2+r3)/2 ha1;
ht1=1/(2(1+r^2-2r w)) ((mb-mc)/(mB-mDStar) (1-r)^2 (w+1)ha1-(mb+mc)/(mB+mDStar) (1+r)^2 (w-1)hv);
ht2=((1-r^2)(w+1))/(2(1+r^2-2r w)) ((mb-mc)/(mB-mDStar) ha1-(mb+mc)/(mB+mDStar) hv);
ht3=-(1/(2(1+r)(1+r^2-2r w)))(2 (mb-mc)/(mB-mDStar) r(w+1)ha1-(mb-mc)/(mB-mDStar) (1+r^2-2r w)(ha3-r ha2)-(mb+mc)/(mB+mDStar) (1+r)^2 hv);

v=(mB+mDStar)/(2Sqrt[mB mDStar]) hv;
a1=((mB+mDStar)^2-qs)/(2Sqrt[mB mDStar](mB+mDStar)) ha1;
a2=(mB+mDStar)/(2Sqrt[mB mDStar]) (ha3+mDStar/mB ha2);
a0=1/(2Sqrt[mB mDStar]) (((mB+mDStar)^2-qs)/(2mDStar) ha1-(mB^2-mDStar^2+qs)/(2mB) ha2-(mB^2-mDStar^2-qs)/(2mDStar) ha3);
t1=1/(2Sqrt[mB mDStar]) ((mB+mDStar)ht1-(mB-mDStar)ht2);
t2=1/(2Sqrt[mB mDStar]) (((mB+mDStar)^2-qs)/(mB+mDStar) ht1-((mB-mDStar)^2-qs)/(mB-mDStar) ht2);
t3=1/(2Sqrt[mB mDStar]) ((mB-mDStar)ht1-(mB+mDStar)ht2-(2(mB^2-mDStar^2))/mB ht3);

alDstarMinus=(8mB^2 magpDstar^2)/(mB+mDStar)^2 Abs[cvl]^2 v^2+((mB+mDStar)^2 (8mDStar^2 qs + \[Lambda]))/(2mDStar^2 qs) Abs[cal]^2 a1^2+(8mB^4 magpDstar^4)/(mDStar^2 (mB+mDStar)^2 qs) Abs[cal]^2 a2^2-(4magpDstar^2 mB^2 (mB^2-mDStar^2-qs))/(mDStar^2 qs) Abs[cal]^2 a1 a2+mLeps[[l]]((32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1+(8(mB+mDStar)(2mDStar^2 (mB^2-mDStar^2)+mB^2 magpDstar^2))/(qs mDStar^2) Re[cal Conjugate[ctl]]a1 t2-(8mB^2 (mB^2-mDStar^2-qs)magpDstar^2)/(qs (mB-mDStar)mDStar^2) Re[cal Conjugate[ctl]]a1 t3-(8mB^2 (mB^2+3mDStar^2-qs)magpDstar^2)/(qs(mB+mDStar)mDStar^2) Re[cal Conjugate[ctl]]a2 t2+(32mB^4 magpDstar^4)/(qs mDStar^2 (mB+mDStar)(mB^2-mDStar^2)) Re[cal Conjugate[ctl]]a2 t3)+mLeps[[l]]^2 ((32mB^2 magpDstar^2)/qs^2 Abs[ctl]^2 t1^2+(2(8mDStar^2 (2(mB^2+mDStar^2)-qs)qs+(4mDStar^2+qs)\[Lambda]))/(qs^2 mDStar^2) Abs[ctl]^2 t2^2+(32mB^4 magpDstar^4)/(qs mDStar^2 (mB^2-mDStar^2)^2) Abs[ctl]^2 t3^2-(16mB^2 magpDstar^2 (mB^2+3mDStar^2-qs))/(qs mDStar^2 (mB^2-mDStar^2)) Abs[ctl]^2 t2 t3);

clDstarMinus=(8magpDstar^2 mB^2)/(mB+mDStar)^2 Abs[cvl]^2 v^2-((mB+mDStar)^2 \[Lambda])/(2mDStar^2 qs) Abs[cal]^2 a1^2-(8 magpDstar^4 mB^4)/((mB+mDStar)^2 mDStar^2 qs) Abs[cal]^2 a2^2+(4magpDstar^2 mB^2 (mB^2-mDStar^2-qs))/(mDStar^2 qs) Abs[cal]^2 a1 a2+mLeps[[l]]((32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1-(8mB^2 (mB+mDStar)magpDstar^2)/(qs mDStar^2) Re[cal Conjugate[ctl]]a1 t2+(8mB^2 (mB^2-mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB-mDStar)) Re[cal Conjugate[ctl]]a1 t3 + (8mB^2 (mB^2+3mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB+mDStar)) Re[cal Conjugate[ctl]]a2 t2-(32mB^4 magpDstar^4)/(qs mDStar^2 (mB+mDStar)(mB^2-mDStar^2)) Re[cal Conjugate[ctl]]a2 t3)+mLeps[[l]]^2 ((32mB^2 magpDstar^2)/qs^2 Abs[ctl]^2 t1^2+(2(4mDStar^2-qs)\[Lambda])/(mDStar^2 qs^2) Abs[ctl]^2 t2^2-(32mB^4 magpDstar^4)/(qs mDStar^2 (mB^2-mDStar^2)^2) Abs[ctl]^2 t3^2+(16mB^2 magpDstar^2 (mB^2+3mDStar^2-qs))/(qs mDStar^2 (mB^2-mDStar^2)) Abs[ctl]^2 t2 t3);

nPrefactor=(\[Tau]B Abs[GFVcb]^2 qs)/(256 \[Pi]^3 mB^2) (1-mLeps[[l]]^2/qs)^2;
branchingRatioMinus=nPrefactor magpDstar (2 alDstarMinus + 2/3 clDstarMinus);

branchingRatioMinus
]
];

bardhand\[CapitalGamma]BDstarl\[Nu]dqsPlus=Function[{cvl,cal,cpl,ctl,qs,l},
Block[
{magpDstar,\[Lambda],
r,w,z,
r11,r21,\[Rho]DstarSq,ha11,
ha1,r1,r2,r3,
hv,ha2,ha3,ht1,ht2,ht3,
v,a0,a1,a2,t1,t2,t3,
alDstarMinus,clDstarMinus,alDstarPlus,clDstarPlus,
nPrefactor,branchingRatioMinus,branchingRatioPlus},

magpDstar=Sqrt[lambda[mB^2,mDStar^2,qs]]/(2mB);
(* Not sure if this is correct *)
\[Lambda]=lambda[mB^2,mDStar^2,qs];

r=mDStar/mB;
w=(mB^2+mDStar^2-qs)/(2mB mDStar);
z=(Sqrt[w+1]-Sqrt[2])/(Sqrt[w+1]+Sqrt[2]);

r11=1.406;
r21=0.853;
\[Rho]DstarSq=1.207;
ha11=0.906;

ha1=ha11 (1-8 \[Rho]DstarSq z + (53\[Rho]DstarSq - 15)z^2-(231\[Rho]DstarSq - 91)z^3);
r1=r11-0.12(w-1)+0.05(w-1)^2;
r2=r21+0.11(w-1)-0.06(w-1)^2;
r3=1.22 -0.052(w-1)+0.026(w-1)^2;

hv=r1 ha1;
ha2=(r2-r3)/(2r) ha1;
ha3=(r2+r3)/2 ha1;
ht1=1/(2(1+r^2-2r w)) ((mb-mc)/(mB-mDStar) (1-r)^2 (w+1)ha1-(mb+mc)/(mB+mDStar) (1+r)^2 (w-1)hv);
ht2=((1-r^2)(w+1))/(2(1+r^2-2r w)) ((mb-mc)/(mB-mDStar) ha1-(mb+mc)/(mB+mDStar) hv);
ht3=-(1/(2(1+r)(1+r^2-2r w)))(2 (mb-mc)/(mB-mDStar) r(w+1)ha1-(mb-mc)/(mB-mDStar) (1+r^2-2r w)(ha3-r ha2)-(mb+mc)/(mB+mDStar) (1+r)^2 hv);

v=(mB+mDStar)/(2Sqrt[mB mDStar]) hv;
a1=((mB+mDStar)^2-qs)/(2Sqrt[mB mDStar](mB+mDStar)) ha1;
a2=(mB+mDStar)/(2Sqrt[mB mDStar]) (ha3+mDStar/mB ha2);
a0=1/(2Sqrt[mB mDStar]) (((mB+mDStar)^2-qs)/(2mDStar) ha1-(mB^2-mDStar^2+qs)/(2mB) ha2-(mB^2-mDStar^2-qs)/(2mDStar) ha3);
t1=1/(2Sqrt[mB mDStar]) ((mB+mDStar)ht1-(mB-mDStar)ht2);
t2=1/(2Sqrt[mB mDStar]) (((mB+mDStar)^2-qs)/(mB+mDStar) ht1-((mB-mDStar)^2-qs)/(mB-mDStar) ht2);
t3=1/(2Sqrt[mB mDStar]) ((mB-mDStar)ht1-(mB+mDStar)ht2-(2(mB^2-mDStar^2))/mB ht3);

alDstarPlus=(8magpDstar^2 mB^2)/(mb+mc)^2 Abs[cpl]^2 a0^2+(32mB^2 magpDstar^2)/qs Abs[ctl]^2 t1^2+(8(mB^2-mDStar^2)^2)/qs Abs[ctl]^2 t2^2-mLeps[[l]]((16 magpDstar^2 mB^2)/((mb+mc)qs) Re[cal Conjugate[cpl]]a0^2-(32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1-(8(mB+mDStar)(mB^2-mDStar^2))/qs Re[cal Conjugate[ctl]]a1 t2)+mLeps[[l]]^2 ((8magpDstar^2 mB^2)/qs^2 Abs[cal]^2 a0^2+(8magpDstar^2 mB^2)/((mB+mDStar)^2 qs) Abs[cvl]^2 v^2+(2(mB+mDStar)^2)/qs Abs[cal]^2 a1^2);

clDstarPlus=-((32mB^2 magpDstar^2)/qs) Abs[ctl]^2 t1^2-(2(4mDStar^2-qs)\[Lambda])/(mDStar^2 qs) Abs[ctl]^2 t2^2+(32mB^4 magpDstar^4)/(mDStar^2 (mB^2-mDStar^2)^2) Abs[ctl]^2 t3^2-(16mB^2 magpDstar^2 (mB^2+3mDStar^2-qs))/(mDStar^2 (mB^2-mDStar^2)) Abs[ctl]^2 t2 t3-mLeps[[l]]((32mB^2 magpDstar^2)/(qs(mB+mDStar)) Re[cvl Conjugate[ctl]]v t1 - (8mB^2 (mB+mDStar)magpDstar^2)/(qs mDStar^2) Re[cal Conjugate[ctl]]a1 t2+(8mB^2 (mB^2-mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB-mDStar)) Re[cal Conjugate[ctl]]a1 t3+(8mB^2 (mB^2+3mDStar^2-qs)magpDstar^2)/(qs mDStar^2 (mB+mDStar)) Re[cal Conjugate[ctl]]a2 t2-(32mB^4 magpDstar^4)/(qs mDStar^2 (mB+mDStar)(mB^2-mDStar^2)) Re[cal Conjugate[ctl]]a2 t3)+mLeps[[l]]^2 (-((8magpDstar^2 mB^2)/((mB+mDStar)^2 qs)) Abs[cvl]^2 v^2+((mB+mDStar)^2 \[Lambda])/(2mDStar^2 qs^2) Abs[cal]^2 a1^2+(8magpDstar^4 mB^4)/(mDStar^2 (mB+mDStar)^2 qs^2) Abs[cal]^2 a2^2-(4magpDstar^2 mB^2)/(mDStar^2 qs^2) (mB^2-mDStar^2-qs)Abs[cal]^2 a1 a2);

nPrefactor=(\[Tau]B Abs[GFVcb]^2 qs)/(256 \[Pi]^3 mB^2) (1-mLeps[[l]]^2/qs)^2;
branchingRatioPlus=nPrefactor magpDstar (2 alDstarPlus + 2/3 clDstarPlus);
branchingRatioPlus
]
];

rDStar=Compile[{{m,_Real},{cV,_Real},{cS,_Real},{cT,_Real}},
Block[{},

NIntegrate[bardhand\[CapitalGamma]BDstarl\[Nu]dqs[1+cV,-(cV+1),-runningS[m] cS,2 runningT[m] cT,qs,3],{qs,mLeps[[3]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]6,*)WorkingPrecision->7,AccuracyGoal->5]/(NIntegrate[bardhand\[CapitalGamma]BDstarl\[Nu]dqs[1,-1,0,0,qs,2],{qs,mLeps[[2]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]6,*)WorkingPrecision->7,AccuracyGoal->5])
],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];

rDStar\[Mu]mLQ=Compile[{{cV,_Real},{cS,_Real},{cT,_Real}},
Block[{runningS,runningT},
runningS=1;
runningT=1;
NIntegrate[bardhand\[CapitalGamma]BDstarl\[Nu]dqs[1+cV,-(cV+1),-runningS cS,2 runningT cT,qs,3],{qs,mLeps[[3]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]6,*)WorkingPrecision->7,AccuracyGoal->5]/(NIntegrate[bardhand\[CapitalGamma]BDstarl\[Nu]dqs[1,-1,0,0,qs,2],{qs,mLeps[[2]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},(*MinRecursion\[Rule]6,*)WorkingPrecision->7,AccuracyGoal->5])
],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];

rDStarLQ=Function[{m,x,y},
Block[
{z,
cV\[Tau],cS\[Tau],cT\[Tau],
cV\[Mu],cS\[Mu],cT\[Mu],
\[CapitalGamma]\[Tau],\[CapitalGamma]\[Mu]},

z=x.ConjugateTranspose[CKM];

cV\[Tau]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[3,2]]x[[;;,3]];
cS\[Tau]=1/(2Sqrt[2]GFVcb) y[[3,2]]/(2m^2) x[[;;,3]];
cT\[Tau]=(-1/4) cS\[Tau];

cV\[Mu]=1/(2Sqrt[2]GFVcb) 1/(2m^2) z[[2,2]]x[[;;,3]];
cS\[Mu]=1/(2Sqrt[2]GFVcb) y[[2,2]]/(2m^2) x[[;;,3]];
cT\[Mu]=(-1/4) cS\[Mu];

\[CapitalGamma]\[Tau]=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDstarl\[Nu]dqs[KroneckerDelta[3,\[Nu]l]+cV\[Tau][[\[Nu]l]],-(cV\[Tau][[\[Nu]l]]+KroneckerDelta[3,\[Nu]l]),-runningS[m] cS\[Tau][[\[Nu]l]],2 runningT[m] cT\[Tau][[\[Nu]l]],qs,3],
{\[Nu]l,1,3}],
{qs,mLeps[[3]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},MinRecursion->9,WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]\[Mu]=NIntegrate[
Sum[
bardhand\[CapitalGamma]BDstarl\[Nu]dqs[KroneckerDelta[2,\[Nu]l]+cV\[Mu][[\[Nu]l]],-(cV\[Mu][[\[Nu]l]]+KroneckerDelta[2,\[Nu]l]),-runningS[m] cS\[Mu][[\[Nu]l]],2 runningT[m] cT\[Mu][[\[Nu]l]],qs,2],
{\[Nu]l,1,3}],
{qs,mLeps[[2]]^2,(mB-mDStar)^2},Method->{Automatic,"SymbolicProcessing"->0},MinRecursion->9,WorkingPrecision->7,AccuracyGoal->5];

\[CapitalGamma]\[Tau]/\[CapitalGamma]\[Mu]]];


CLL[m_,x_]:=With[{z=x.ConjugateTranspose[CKM],mt=162.3,\[Alpha]=1/127},mt^2/(8\[Pi] \[Alpha] (m)^2) Abs[z[[2,3]]]^2-1/(64\[Pi] \[Alpha]) Sqrt[2]/(GF (m)^2) 1/(CKM[[3,3]]Conjugate[CKM[[3,2]]]) (Sum[x[[j,3]]Conjugate[x[[j,2]]],{j,3}])(Abs[z[[2,1]]]^2+Abs[z[[2,2]]]^2+Abs[z[[2,3]]]^2)]//Re
CLR[m_,x_,y_]:=With[{z=x.ConjugateTranspose[CKM],mt=162.3,\[Alpha]=1/127},
mt^2/(16\[Pi] \[Alpha] (m)^2) Abs[y[[2,3]]]^2-1/(64\[Pi] \[Alpha]) Sqrt[2]/(GF (m)^2) 1/(CKM[[3,3]]Conjugate[CKM[[3,2]]]) (Sum[x[[j,3]]Conjugate[x[[j,2]]],{j,3}])(Abs[y[[2,1]]]^2+Abs[y[[2,2]]]^2+Abs[y[[2,3]]]^2)]//Re
hh[qs_,m_]:=With[{x=4m^2/qs},\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{
RowBox[{
FractionBox[
RowBox[{"-", "4"}], "9"], 
RowBox[{"(", 
RowBox[{
RowBox[{"Log", "[", 
FractionBox[
SuperscriptBox["m", "2"], "qs"], "]"}], "-", 
FractionBox["2", "3"], "-", "x"}], ")"}]}], "-", 
RowBox[{
FractionBox["4", "9"], 
RowBox[{"(", 
RowBox[{"2", "+", "x"}], ")"}], " ", 
SqrtBox[
RowBox[{"x", "-", "1"}]], 
RowBox[{"ArcTan", "[", 
FractionBox["1", 
SqrtBox[
RowBox[{"x", "-", "1"}]]], "]"}]}]}], 
RowBox[{"x", ">", "1"}]},
{
RowBox[{
RowBox[{
FractionBox[
RowBox[{"-", "4"}], "9"], 
RowBox[{"(", 
RowBox[{
RowBox[{"Log", "[", 
FractionBox[
SuperscriptBox["m", "2"], "qs"], "]"}], "-", 
FractionBox["2", "3"], "-", "x"}], ")"}]}], "-", 
RowBox[{
FractionBox["4", "9"], 
RowBox[{"(", 
RowBox[{"2", "+", "x"}], ")"}], 
SqrtBox[
RowBox[{"1", "-", "x"}]], 
RowBox[{"(", 
RowBox[{"Log", "[", 
RowBox[{
FractionBox[
RowBox[{"1", "+", 
SqrtBox[
RowBox[{"1", "-", "x"}]]}], 
SqrtBox["x"]], "-", 
FractionBox[
RowBox[{"I", " ", "\[Pi]"}], "2"]}], "]"}], ")"}]}]}], 
RowBox[{"x", "<=", "1"}]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False]\)];
hh[qs_,0]:=8/27; (* Took limit of above *)
(* Many values for constants taken from table V in 1306.2384 *)
dBrdqsRK[m_,x_,y_,qs_,l_]:=Block[{c9,c10,nksq,\[Tau]Bd=1.519 10^-12,\[Alpha]=1/128.957,\[Lambda]q,\[Lambda]B,mK=0.4976,\[Phi]9,\[Phi]10,mB=5.279,mb=4.18,mc=1.275,
a0,aPlus,t0,mPlus,tPlus,zz,fPlus,f0,c9sm,c10sm,Y,
c1=-0.257,c2=1.009,c3=-0.0050,c4=-0.078,c5=0, c6=0.001},

(* From Becerivic et al. [1608.07583] *)
Y=4/3 c3+64/9 c5+64/27 c6-1/2 hh[qs,0](c3+4/3 c4+16c5+64/3 c6)+hh[qs,mc](4/3 c1 + c2 + 6c3 + 60c5) - 1/2 hh[qs ,mb] (7c3 + 4/3 c4 + 76 c5 +64/3 c6);
c9sm=4.211+Y;
c10sm=-4.103;

c9=1/2 (CLL[m,x]+CLR[m,x,y])KroneckerDelta[l,2]+c9sm;
c10=1/2 (-CLL[m,x]+CLR[m,x,y])KroneckerDelta[l,2]+c10sm;
\[Lambda]q=(qs-4(mLeps[[l]])^2)(qs);
\[Lambda]B=(qs-(mB+mK)^2)(qs-(mB-mK)^2);

(* From 1411.3161 *)
Subscript[a0, 0]=0.54; Subscript[a0, 1]=-1.91; Subscript[a0, 2]=1.83; Subscript[a0, 3]=-0.02;
Subscript[aPlus, 0]=0.43; Subscript[aPlus, 1]=-0.67;Subscript[aPlus, 2]=-1.12;

(* From 1306.2384 *)
t0=(mB+mK)(Sqrt[mB]-Sqrt[mK])^2;
mPlus=mB+0.046;
tPlus=(mB+mK)^2;
zz=(Sqrt[tPlus-qs]-Sqrt[tPlus-t0])/(Sqrt[tPlus-qs]+Sqrt[tPlus-t0]);
fPlus=1/(1-qs/mPlus^2) (Subscript[aPlus, 0]+Subscript[aPlus, 1] zz+Subscript[aPlus, 2] zz^2+zz^3/3 (-Subscript[aPlus, 1]+2Subscript[aPlus, 2]));
f0=Sum[Subscript[a0, j] zz^j,{j,0,3}];

nksq=(\[Tau]Bd/hbar) (\[Alpha]^2 GF^2 Abs[CKM[[3,3]]Conjugate[CKM[[3,2]]]]^2)/(512\[Pi]^5 mB^3) Sqrt[\[Lambda]q \[Lambda]B]/qs;\[Phi]9=1/2 Abs[fPlus]^2 \[Lambda]B(1-\[Lambda]q/(3qs^2));\[Phi]10=1/2 Abs[f0]^2 4mLeps[[l]]^2 (mB^2-mK^2)^2/qs+1/2 Abs[fPlus]^2 \[Lambda]B(1-(4mLeps[[l]]^2)/qs-\[Lambda]q/(3qs^2));Abs[nksq](\[Phi]9 Abs[c9]^2+\[Phi]10 Abs[c10]^2)];
(*rKLQ=Compile[{{m,_Real},{x,_Real,2},{y, _Real,2}},
With[{mK=0.4976,mB=5.279},
NIntegrate[dBrdqsRK[m,x,y,qs,2],{qs,4m\[Mu]^2,(mB-mK)^2},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]/NIntegrate[dBrdqsRK[m,x,y,qs,1],{qs,4me^2,(mB-mK)^2},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]
],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];*)
rKLQ=Compile[{{m,_Real},{x,_Real,2},{y, _Real,2}},
With[{mK=0.4976,mB=5.279},
NIntegrate[dBrdqsRK[m,x,y,qs,2],{qs,1,6},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]/NIntegrate[dBrdqsRK[m,x,y,qs,1],{qs,1,6},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]
],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];


dBrdqsRKOPS[cll_,clr_,qs_,l_]:=Block[{c9,c10,nksq,\[Tau]Bd=1.519 10^-12,\[Alpha]=1/128.957,\[Lambda]q,\[Lambda]B,mK=0.4976,\[Phi]9,\[Phi]10,mB=5.279,mb=4.18,mc=1.275,
a0,aPlus,t0,mPlus,tPlus,zz,fPlus,f0,c9sm,c10sm,Y,
c1=-0.257,c2=1.009,c3=-0.0050,c4=-0.078,c5=0, c6=0.001},

(* From Becerivic et al. [1608.07583] *)
Y=4/3 c3+64/9 c5+64/27 c6-1/2 hh[qs,0](c3+4/3 c4+16c5+64/3 c6)+hh[qs,mc](4/3 c1 + c2 + 6c3 + 60c5) - 1/2 hh[qs ,mb] (7c3 + 4/3 c4 + 76 c5 +64/3 c6);
c9sm=4.211+Y;
c10sm=-4.103;

c9=1/2 (cll+clr)KroneckerDelta[l,2]+c9sm;
c10=1/2 (-cll+clr)KroneckerDelta[l,2]+c10sm;
\[Lambda]q=(qs-4(mLeps[[l]])^2)(qs);
\[Lambda]B=(qs-(mB+mK)^2)(qs-(mB-mK)^2);

(* From 1411.3161 *)
Subscript[a0, 0]=0.54; Subscript[a0, 1]=-1.91; Subscript[a0, 2]=1.83; Subscript[a0, 3]=-0.02;
Subscript[aPlus, 0]=0.43; Subscript[aPlus, 1]=-0.67;Subscript[aPlus, 2]=-1.12;

(* From 1306.2384 *)
t0=(mB+mK)(Sqrt[mB]-Sqrt[mK])^2;
mPlus=mB+0.046;
tPlus=(mB+mK)^2;
zz=(Sqrt[tPlus-qs]-Sqrt[tPlus-t0])/(Sqrt[tPlus-qs]+Sqrt[tPlus-t0]);
fPlus=1/(1-qs/mPlus^2) (Subscript[aPlus, 0]+Subscript[aPlus, 1] zz+Subscript[aPlus, 2] zz^2+zz^3/3 (-Subscript[aPlus, 1]+2Subscript[aPlus, 2]));
f0=Sum[Subscript[a0, j] zz^j,{j,0,3}];

nksq=(\[Tau]Bd/hbar) (\[Alpha]^2 GF^2 Abs[CKM[[3,3]]Conjugate[CKM[[3,2]]]]^2)/(512\[Pi]^5 mB^3) Sqrt[\[Lambda]q \[Lambda]B]/qs;\[Phi]9=1/2 Abs[fPlus]^2 \[Lambda]B(1-\[Lambda]q/(3qs^2));\[Phi]10=1/2 Abs[f0]^2 4mLeps[[l]]^2 (mB^2-mK^2)^2/qs+1/2 Abs[fPlus]^2 \[Lambda]B(1-(4mLeps[[l]]^2)/qs-\[Lambda]q/(3qs^2));Abs[nksq](\[Phi]9 Abs[c9]^2+\[Phi]10 Abs[c10]^2)];

rK=Compile[{{cll,_Real},{clr,_Real}},
With[{mK=0.4976,mB=5.279},
NIntegrate[dBrdqsRKOPS[cll,clr,qs,2],{qs,1,6},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]/NIntegrate[dBrdqsRKOPS[cll,clr,qs,1],{qs,1,6},Method->{Automatic,"SymbolicProcessing"->0},WorkingPrecision->7,AccuracyGoal->5]
],
CompilationTarget->"C",
RuntimeAttributes->{Listable},
Parallelization->True
];


(*need to add condition*)
\[Delta]a\[Mu][m_,x_,y_]:=With[{z=x.ConjugateTranspose[CKM], mu={0.0022,1.25,160}},
Sum[(m\[Mu] mu[[q]])/(4\[Pi]^2 m^2) (Log[m^2/mu[[q]]^2]-7/4)Re[-Conjugate[y[[2,q]]] z[[2,q]]],{q,1,3}]-m\[Mu]^2/(32\[Pi]^2 m^2) (Sum[Abs[z[[2,i]]]^2,{i,3}]+Sum[Abs[y[[2,j]]]^2,{j,3}])]


<<HypothesisTesting`
confidenceLevel2d[x_,\[Mu]x_,\[Sigma]x_,y_,\[Mu]y_,\[Sigma]y_]:=Block[{\[Chi]2},
\[Chi]2=((x-\[Mu]x)/\[Sigma]x)^2+((y-\[Mu]y)/\[Sigma]y)^2;
\[Chi]2
(*(1-Min[GammaRegularized[1/2,\[Chi]2/2],GammaRegularized[1/2,0,\[Chi]2/2]])*)]
clrdrdstar[x_,y_]:=confidenceLevel2d[x,0.397,Sqrt[0.04^2+0.028^2],y,0.311,0.016]//N;
