(* ::Package:: *)

Print["PhysicsWorkbench Version: 1.0 (8 Feb 2017)."];
Print["Author: John Gargalionis"];
(*
Print["Please cite: arXiv:06051991"];
Print["The package Physics Workbench is written for Mathematica 10 and higher."];
Print["Physics Workbench provides various useful constants, matrices and functions used in particle phenomenology. All masses are given in GeV. Quark masses are purposefully omitted to draw attention to running effects."];
*)


BeginPackage["PhysicsWorkbench`"];

{CKM, \[Theta]w, sw2, mLeps, me, m\[Mu], m\[Tau], mW, mZ, mB, mD, mDStar, \[Tau]B, vev, hbar, GF, GFVcb}

EndPackage[];


(* Definition of CKM matrix *)
CKM=With[{\[Lambda]=0.2257,A=0.814,\[Rho]=0.135,\[Eta]=0.349},
({
 {1-\[Lambda]^2/2, \[Lambda], A \[Lambda]^3 (\[Rho]-I \[Eta])},
 {-\[Lambda], 1-\[Lambda]^2/2, A \[Lambda]^2},
 {A \[Lambda]^3 (1-\[Rho]-I \[Eta]), -A \[Lambda]^2, 1}
})];
CKM//MatrixForm

(* Weak mixing angle *)
\[Theta]w = 28.76 Degree;
sw2 = 0.23142; (* get value from PDG *)

(* Charged lepton masses *)
mLeps={5.11 10^-4,0.106,1.777};
me = mLeps[[1]];
m\[Mu] = mLeps[[2]];
m\[Tau] = mLeps[[3]];

(* quark masses *)
md={0.005,0.095,4.8};  (* index: quark *)
mu={0.0022,1.25,174.2}; (* index quark *)
mb=4.2;
mc=1.25;
mt=174.2;

(* Gauge boson masses *)
mW = 80.403;
mZ = 91.188;

(* Meson properties *)
mB = 5.26; 
mD = 1.86;
mDStar = 2.01;
mK = 493.677 10^-3;
\[Tau]B = (1.638 10^-12)/(6.58 10^-25); (* lifetime of B meson in GeV^-1 *)

(* Other constants *)
vev = 246; (* Higgs VEV in GeV *)
hbar = 6.582119*10^-25; (* GeV sec *)
GF = 1.16638*10^-5;(* Inverse GeV squared *)
GFVcb = (1.16638 10^-5) (41.1 10^-3);
