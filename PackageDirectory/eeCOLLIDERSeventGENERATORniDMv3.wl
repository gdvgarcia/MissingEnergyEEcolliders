(* ::Package:: *)

BeginPackage["eeCOLLIDERSeventGENERATORniDMv3`"];


(* ::Title::Closed:: *)
(*Public*)


(* ::Subsubsection::Closed:: *)
(*Introduction text*)


Print["
****************************************  Electron-positron Colliders Event Generator for niDM Package  *****************************************  

\nThis package allows the user to generate a variety of high energy events at electron-positron colliders in the context of the Standard Model supplemen-
\nted by not-so-inelastic Dark Matter model (DM) with a dark photon mediator. DM notation will refer to the lighest dark fermionic state \\chi- (\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)), e.g. 
\nmDM = \!\(\*SubscriptBox[\(m\), \(\[Chi]\)]\) = mass of chi-. Zp (Z') will refer to the dark photon.
\n
\n\tCurrently available events: {\"e- e+ -> \[Gamma] Z'\"}. ***We will be using the place holder EVENT_TYPE \[Element] {\"e- e+ -> \[Gamma] Z'\"}
\n
\n\tFunctions available (angles in degrees and energies/masses/momenta in GeV - velocities with no units - v = \!\(\*SubscriptBox[\(v\), \(km\)]\)/\!\(\*SubscriptBox[\(c\), \(light\)]\) ): 

\t\t\tFullEventGenerator[EVENT_TYPE,ModelParameters,\[Theta]\[Gamma]minCM,\[Theta]\[Gamma]maxCM,EnergyInCM,VelocityOfCM,n,nPhotons,nChis,verbose]

\t\t\tOnlySMEventSelected[EVENT_TYPE,EventsWithFullInfo,verbose]

\t\t\tFullEventsDisplay[EVENT_TYPE,FullEvent,AxisBool,HidePhotonFlag]

\t\t\tSMEventsDisplay[EVENT_TYPE,SMEvent,AxisBool]

\t\t\t\[Sigma]WithinAngularInterval[EVENT_TYPE,\!\(\*SuperscriptBox[\(EnergyInCM\), \(2\)]\),\[Theta]\[Gamma]minCM,\[Theta]\[Gamma]maxCM,\[Epsilon],mZp]

\n
\nCheck function::usage for more information on each function.
\n
\n
\n							*** Attention: Do not define variables with the same name as any of the functions above: FullEventGenerator, OnlySMEventSelected, etc! ***
\n																					                                                      By Giovani Dalla Valle Garcia (KIT - DAAD), Karlsruhe 2024\n\n"]


(* ::Subsubsection::Closed:: *)
(*Functions usage *)


FullEventGenerator::usage = " returns all information about n events generated.\n\n\tCall it with as follows:\n\t\tFullEventGenerator[EVENT_TYPE(*in {\"e- e+ -> \[Gamma] Z'\"}),ModelParameters(*{mDM,mZp,\[Alpha]D,\[Epsilon],\[Delta]y,\[CapitalDelta]f}*),\[Theta]\[Gamma]minCM(*minimal photon CMangle with the Z-axis in degrees*),\[Theta]\[Gamma]maxCM(*maximal photon CMangle with the Z-axis in degrees*),EnergyInCM(*CM energy in GeV*),VelocityOfCM(*velocity of the CM normalized by the speed of light*),n,nPhotons(*number of points to scan the collision distribution*),nChis(*number of points to scan the chi 3-body decay distribution*),verbose]";


OnlySMEventSelected::usage = " selects only the SM information from the list of full events information EventsWithFullInfo.\n\n\tCall it with as follows:\n\t\tOnlySMEventSelected[EVENT_TYPE(*in {\"e- e+ -> \[Gamma] Z'\"}),EventsWithFullInfo,verbose]";


FullEventsDisplay::usage = " displays a full event graphically - (assuming the second daugther of the dark photon is a chi+ in order for to decay it).\n\n\tCall it with as follows:\n\t\tFullEventsDisplay[EVENT_TYPE(*in {\"e- e+ -> \[Gamma] Z'\"}),FullEvent(*event in the format of the output from tFullEventGenerator[][[i]]*),AxisBool(*0 no axis displayed, 1 axis displayed*),HidePhotonFlag(*0 no photon displayed, 1 photon displayed*)]";


SMEventsDisplay::usage = "  displays the SM part of an event graphically - (assuming the second daugther of the dark photon is a chi+ in order for to decay it).\n\n\tCall it with as follows:\n\t\tSMEventsDisplay[EVENT_TYPE(*in {\"e- e+ -> \[Gamma] Z'\"}),SMEvent(*event in the format of the output from OnlySMEventSelected[][[i]]*),AxisBool(*0 no axis displayed, 1 axis displayed*)]";


\[Sigma]WithinAngularInterval::usage = " returns the cross-section for producing a dark photon together with a photon between the angles \[Theta]\[Gamma]minCM and \[Theta]\[Gamma]maxCM.\n\n\tCall it with as follows:\n\t\t\[Sigma]WithinInterval[\"e- e+ -> \[Gamma] Z'\",\!\(\*SuperscriptBox[\(EnergyInCM\), \(2\)]\)(*in \!\(\*SuperscriptBox[\(GeV\), \(2\)]\)*),\[Theta]\[Gamma]minCM(*in degrees*),\[Theta]\[Gamma]maxCM(*in degrees*),\[Epsilon],mZp]";


(* ::Text:: *)
(*It also defines which functions are public.*)


(* ::Title::Closed:: *)
(*Private*)


Begin["`Private`"];


removeQ=DialogInput[Column[{"Final states of the simulation should...",Row[{Button["include pions.",pionsFlag=True;DialogReturn["Yes"]],Button["not include pions.",pionsFlag=False;DialogReturn["Yes"]]}]}]];


If[!pionsFlag,Print["Procedding with no pions...\n"];,Print["Procedding with pions...\n"];];


removeQ=DialogInput[Column[{"The simulation should be...",Row[{Button["compiled in C.",compilationFlag=True;DialogReturn["Yes"]],Button["not be compiled in C(slow simulations).",compilationFlag=False;DialogReturn["Yes"]]}]}]]


If[!compilationFlag,0;,Needs["CCompilerDriver`"];<<CompiledFunctionTools`;If[CCompilers[]==={},Print["ERROR in checking for C compiler drivers. Event generation may not work or will be considerable slower. Consider fixing the problem with: Needs[\"CCompilerDriver`\"].\n\n"],PrintTemporary["C compiler drivers already loaded. Compiling functions..."]];];


initialTime=Now;


(* ::Section::Closed:: *)
(*Pre definitions*)


(* ::Subsection::Closed:: *)
(*Defining SM variables*)


ElectronMass=0.5109989461*10^-3(*GeV*);
MuonMass=0.1056583745(*GeV*);
PionMass=0.13957061(*GeV*);
ElectromagneticCoupling=Sqrt[(4\[Pi])/128](*try to give a simple estimate for the elecmag coupling at 10 GeV - see below*);
SpeedOfLight=299792458(*m/s*);


(* ::Subsection::Closed:: *)
(*Formulas definitions*)


(* ::Subsubsection::Closed:: *)
(*Decay width functions*)


m23MIN3bodypDFDM2M[M_,m1_,m2_,m3_,m12_]=1/(2m12) (M^2 (-m1^2+m12+m2^2)-(m12-m2^2) (m12-m3^2)+m1^2 (m12+m3^2)-Sqrt[(m1^2-m12)^2-2 (m1^2+m12) m2^2+m2^4] Sqrt[(M^2-m12)^2-2 (M^2+m12) m3^2+m3^4]);
m23MAX3bodypDFDM2M[M_,m1_,m2_,m3_,m12_]=1/(2m12) (M^2 (-m1^2+m12+m2^2)-(m12-m2^2) (m12-m3^2)+m1^2 (m12+m3^2)+Sqrt[(m1^2-m12)^2-2 (m1^2+m12) m2^2+m2^4] Sqrt[(M^2-m12)^2-2 (M^2+m12) m3^2+m3^4]);

\[CapitalGamma]3bodypDFDM2M[process_,M_,m1_,m2_,m3_,mZp_,\[Alpha]D_,\[Epsilon]_,\[Delta]y_,\[CapitalDelta]f_]:=If[M>m1+m2+m3,1/(2\[Pi])^3 1/(32M^3)*Sfactor[process]*NIntegrate[Msquared3bodypDFDM2M[process,M,m1,m2,m3,mZp,\[Alpha]D,\[Epsilon],\[Delta]y,\[CapitalDelta]f,m12,m23],{m12,(m1+m2)^2,(M-m3)^2},{m23,m23MIN3bodypDFDM2M[M,m1,m2,m3,m12],m23MAX3bodypDFDM2M[M,m1,m2,m3,m12]},Method->(*intmethod[process]*)"AdaptiveQuasiMonteCarlo"],0];

\[CapitalGamma]2bodypDFDM2M[process_,M_,m1_,m2_,mZp_,\[Alpha]D_,\[Epsilon]_,\[Delta]_,\[CapitalDelta]_]:=If[m1+m2>=M,0,1/(16\[Pi] ) Sqrt[(M^2-(m1+m2)^2)(M^2-(m1-m2)^2)]/M^3 Sfactor[process]*Msquared2bodypDFDM2M[process,M,m1,m2,mZp,\[Alpha]D,\[Epsilon],\[Delta],\[CapitalDelta]]//Simplify];


(* ::Subsubsection::Closed:: *)
(*e- e+ -> Z' \[Gamma]*)


CrossSectionElectronPositronCollisionWithinAngularInterval[s_,thetaCMmin_,thetaCMmax_,eps_,mZp_]:=Module[{tempFunc},
tempFunc[x_]:=ElectronPositronCollisionDistributionNoEpsCompiled[s,x,mZp];
eps^2 NIntegrate[Hold[tempFunc[x]],{x,(\[Pi]/180)* thetaCMmin,(\[Pi]/180) *thetaCMmax}]
]
\[Sigma]WithinAngularInterval["e- e+ -> \[Gamma] Z'",s_,thetaCMmin_,thetaCMmax_,eps_,mZp_]:=CrossSectionElectronPositronCollisionWithinAngularInterval[s,thetaCMmin,thetaCMmax,eps,mZp]


(* ::Subsubsection::Closed:: *)
(*Z'  decays*)


(* ::Text:: *)
(*Z' -> SubMinus[\[Chi]]  SubMinus[\[Chi]]*)


partlist["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)"]={"Z'","\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)"};
Sfactor["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)"]=1/2;
Msquared2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)",mZp_,mM_,placeHolder_,mZpPlaceHolder_,\[Alpha]D_,\[Epsilon]PlaceHolder_,\[Delta]_,\[CapitalDelta]_]=(16 \[Pi] \[Alpha]D \[Delta]^2 \[CapitalDelta]^2 (mZp^2-4 mM^2))/(3 (\[Delta]+2)^2 (\[CapitalDelta]+2)^2);


(* ::Text:: *)
(*Z' ->  SubMinus[\[Chi]]  SubPlus[\[Chi]]*)


partlist["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"]={"Z'","\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"};
Sfactor["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"]=1;
Msquared2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",mZp_,mM_,mP_,mZpPlaceHolder_,\[Alpha]D_,\[Epsilon]PlaceHolder_,\[Delta]_,\[CapitalDelta]_]=-1/(3 mZp^2) 8 \[Pi] \[Alpha]D (1-(\[Delta]^2 \[CapitalDelta]^2)/((\[Delta]+2)^2 (\[CapitalDelta]+2)^2)) (mM^4+mM^2 (mZp^2-2 mP^2)-6 mM mP mZp^2+mP^4+mP^2 mZp^2-2 mZp^4);


(* ::Text:: *)
(*Z' -> SubPlus[\[Chi]] SubPlus[\[Chi]]*)


partlist["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"]={"Z'","\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"};
Sfactor["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"]=1/2;
Msquared2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",mZp_,mP_,placeHolder_,mZpPlaceHolder_,\[Alpha]D_,\[Epsilon]PlaceHolder_,\[Delta]_,\[CapitalDelta]_]=(16 \[Pi] \[Alpha]D \[Delta]^2 \[CapitalDelta]^2 (mZp^2-4 mP^2))/(3 (\[Delta]+2)^2 (\[CapitalDelta]+2)^2);


(* ::Text:: *)
(*Decay widths*)


(* ::Input:: *)
(*(*List has the from {mDM,mZp,\[Alpha]D,\[Epsilon],\[Theta],\[Delta]y,\[CapitalDelta]f}*)*)


\[CapitalGamma]["Z'",List_]:=\[CapitalGamma]2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)",List[[2]],List[[1]],List[[1]],List[[2]],List[[3]],1,List[[5]],List[[6]]]+\[CapitalGamma]2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List[[2]],List[[1]],(1+List[[6]])*List[[1]],List[[2]],List[[3]],1,List[[5]],List[[6]]]+\[CapitalGamma]2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List[[2]],(1+List[[6]])*List[[1]],(1+List[[6]])*List[[1]],List[[2]],List[[3]],1,List[[5]],List[[6]]];
\[CapitalGamma]["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)",List_]:=\[CapitalGamma]2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)",List[[2]],List[[1]],List[[1]],List[[2]],List[[3]],1,List[[5]],List[[6]]];
\[CapitalGamma]["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List_]:=\[CapitalGamma]2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List[[2]],List[[1]],(1+List[[6]])*List[[1]],List[[2]],List[[3]],1,List[[5]],List[[6]]];
\[CapitalGamma]["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List_]:=\[CapitalGamma]2bodypDFDM2M["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List[[2]],(1+List[[6]])*List[[1]],(1+List[[6]])*List[[1]],List[[2]],List[[3]],1,List[[5]],List[[6]]];


(* ::Subsubsection::Closed:: *)
(*SubPlus[\[Chi]] decays*)


(* ::Text:: *)
(*SubPlus[\[Chi]]-> SubMinus[\[Chi]] SuperPlus[l] SuperMinus[l]*)


partlist["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(+\)]\)"]={"\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)"(*m1*),"\!\(\*SuperscriptBox[\(l\), \(-\)]\)"(*m2*),"\!\(\*SuperscriptBox[\(l\), \(+\)]\)"(*m3*)};
Sfactor["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(+\)]\)"]=1;
Msquared3bodypDFDM2M["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(+\)]\)",M_,m1_,m2_,m3_,mZp_,\[Alpha]D_,\[Epsilon]_,\[Delta]_,\[CapitalDelta]_,m12_,m23_]=1/(mZp^4 (m23-mZp^2)^2) 1.7231563788153554` (1/(mZp^2-8315.178393760001`) 0.12132142105910404` ((2 m1^2+m2^2+m3^2-2 m12-m23-2 m2 m3) M^2-2 m1 (m2^2-4 m3 m2+m3^2-m23) M+2 m12^2-2 m12 m2^2+m23^2+2 m2^2 m3^2-2 m12 m3^2-m23 m3^2+m1^2 (m2^2-2 m3 m2+m3^2-2 m12-m23)-m2^2 m23+2 m12 m23+2 m2 m23 m3) mZp^6-2.866054752661671` ((2 m1^2+m2^2+m3^2-2 m12-m23-2 m2 m3) M^2-2 m1 (m2^2-4 m3 m2+m3^2-m23) M+2 m12^2-2 m12 m2^2+m23^2+2 m2^2 m3^2-2 m12 m3^2-m23 m3^2+m1^2 (m2^2-2 m3 m2+m3^2-2 m12-m23)-m2^2 m23+2 m12 m23+2 m2 m23 m3) mZp^4-4 (M-m1) (-((0.06066071052955202` (m2-m3) (m2 M^3+m1 m2 M^2+(m2^3+m3 m2^2-m23 m2+m1^2 m3-m12 (m2+m3)) M+m1 (m3 (m1^2-m23+m3 (m2+m3))-m12 (m2+m3))) mZp^2)/(mZp^2-8315.178393760001`))+1.4330273763308354` (m2-m3) (m2 M^3+m1 m2 M^2+(m2^3+m3 m2^2-m23 m2+m1^2 m3-m12 (m2+m3)) M+m1 (m3 (m1^2-m23+m3 (m2+m3))-m12 (m2+m3)))+1/(mZp^2-8315.178393760001`)^2 0.30511302550508157` (0.5021039710465376` M m2^4-0.5021039710465376` M (m3^2+m23) m2^2+M^3 (0.5021039710465376` m2+0.49789602895346247` m3) m2+M^2 m1 (0.5021039710465376` m2+0.49789602895346247` m3) m2-0.49789602895346247` M (m1^2+m23) m3 m2-0.5021039710465376` M m1^2 m3^2-0.5021039710465376` M m12 (m2^2-1.` m3^2)-0.5021039710465376` m1 m12 (m2^2-1.` m3^2)-0.5021039710465376` m1 m3 (1.` m3^3-1.` m2^2 m3-1.` m23 m3-0.9916193809734177` m2 m23+m1^2 (0.9916193809734177` m2+1.` m3))) mZp^4) mZp^2+(M-m1)^2 (1/(mZp^2-8315.178393760001`)^2 0.30511302550508157` (-0.5021039710465376` m2^4+0.5021039710465376` (2 m3^2+m23) m2^2+0.995792057906925` m23 m3 m2-0.5021039710465376` m3^2 (m3^2-1.` m23)) mZp^4+(0.06066071052955202` (m2^2+2 m3 m2+m3^2-m23) (m2-m3)^2 mZp^2)/(mZp^2-8315.178393760001`)-1.4330273763308354` (m2-m3)^2 (m2^2+2 m3 m2+m3^2-m23)) (M^2+2 m1 M+m1^2-m23)-1/(mZp^2-8315.178393760001`)^2 0.6127938468964998` ((1.` m1^2+0.5` m2^2+0.5` m3^2-1.` m12-0.4999999999999999` m23+0.9916193809734177` m2 m3) M^2+m1 (-1.` m2^2-3.9664775238936714` m3 m2-1.` m3^2+1.` m23) M+1.` m12^2-0.9999999999999998` m12 m2^2+0.5` m23^2+1.` m2^2 m3^2-0.9999999999999998` m12 m3^2-0.4999999999999999` m23 m3^2+m1^2 (0.5` m2^2+0.9916193809734178` m3 m2+0.5` m3^2-1.` m12-0.5` m23)-0.4999999999999999` m2^2 m23+1.` m12 m23-0.9916193809734177` m2 m23 m3) mZp^8) \[Alpha]D (1-(\[Delta]^2 \[CapitalDelta]^2)/((\[Delta]+2)^2 (\[CapitalDelta]+2)^2)) \[Epsilon]^2;


(* ::Text:: *)
(*SubPlus[\[Chi]]-> SubMinus[\[Chi]] SuperPlus[pion] SuperMinus[pion]*)


Fpi2VDM[m23_]=1.1265731376593378` Abs[((0.6547926192842074` +0.00001179288147985907` I) (1+(0.001004506501192573` -0.000011050017200740626` I)/((0.6125410225` -0.006644698499999999` I)-m23)))/((0.-0.23479121194530173` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+1.8783296955624138` (0.04271767567565097` (0.6010280675999999` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.4990469642661888`))+0.6010280675999999`)+(0.020154984130451314` +0.0169763089414969` I)/((0.-0.04949270045567086` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+0.3959416036453669` (0.040860712191285925` (5.080516` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.8774686004328772`))+5.080516`)+(0.044963937180827813` +0.24598370047761595` I)/((0.-0.17569710967673408` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+1.4055768774138726` (0.041181101512326464` (3.463321` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.8134044037295003`))+3.463321`)-(0.3183349993211615` +0.22649769503854025` I)/((0.-0.3016808572181247` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+2.4134468577449977` (0.041623970301305685` (2.2290490000000003` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.7383136955592318`))+2.2290490000000003`)];
partlist["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)"]={"\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)"(*m1*),"\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)"(*m2*),"\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)"(*m3*)};
Sfactor["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)"]=1;
Msquared3bodypDFDM2M["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)",M_,m1_,m2_,m3_,mZp_,\[Alpha]D_,eps_,\[Delta]_,\[CapitalDelta]_,m12_,m23_]=1/(mZp^4 (m23-mZp^2)^2) 0.6173325661353781` \[Alpha]D eps^2 (1-(\[Delta]^2 \[CapitalDelta]^2)/((\[Delta]+2)^2 (\[CapitalDelta]+2)^2)) Fpi2VDM[m23]^2 (2 mZp^4 (M^4+M^2 (2 m1^2-4 m12-m23)+2 M m1 (2 m2^2-m23+2 m3^2)+m1^4-m1^2 (4 m12+m23)+4 (m12^2-m12 (m2^2-m23+m3^2)+m2^2 m3^2))+2 (M-m1) (m2^2-m3^2) ((M-m1) (m2^2-m3^2) (M^2+2 M m1+m1^2-m23)+2 mZp^2 (M^3+M^2 m1+M (m1^2-2 m12+2 m2^2-m23)+m1 (m1^2-2 m12-m23+2 m3^2))))//Simplify;


(* ::Text:: *)
(*Decays widths*)


(* ::Input:: *)
(*(*List has the from {mDM,mZp,\[Alpha]D,\[Epsilon],\[Delta]y,\[CapitalDelta]f}*)*)


\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(e\), \(-\)]\)\!\(\*SuperscriptBox[\(e\), \(+\)]\)",List_]:=\[CapitalGamma]3bodypDFDM2M["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(+\)]\)",List[[1]]*(1+List[[6]]),List[[1]],ElectronMass,ElectronMass,List[[2]],List[[3]],List[[4]],List[[5]],List[[6]]]
\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Mu]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Mu]\), \(+\)]\)",List_]:=\[CapitalGamma]3bodypDFDM2M["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(-\)]\)\!\(\*SuperscriptBox[\(l\), \(+\)]\)",List[[1]]*(1+List[[6]]),List[[1]],MuonMass,MuonMass,List[[2]],List[[3]],List[[4]],List[[5]],List[[6]]]
\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)",List_]:=\[CapitalGamma]3bodypDFDM2M["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)",List[[1]]*(1+List[[6]]),List[[1]],PionMass,PionMass,List[[2]],List[[3]],List[[4]],List[[5]],List[[6]]]


\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",List_]:=\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(e\), \(-\)]\)\!\(\*SuperscriptBox[\(e\), \(+\)]\)",List]+\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Mu]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Mu]\), \(+\)]\)",List]+\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)",List];


(* ::Subsection::Closed:: *)
(*Event Visualization functions*)


(* ::Subsubsection::Closed:: *)
(*Visualization ThreeBody*)


ThreeBodyEventsDisplay[MotherParticleEvent_,DaughtersEvent_,AxisBool_,HideMotherMomenumFlag_]:=Module[{VertexT,norm,vertex2Visualize,eps,axis,motherVisualize,daugther1Visualize,daugther2Visualize,daugther3Visualize},

VertexT=DaughtersEvent[[2]];
norm=0.3*Norm[MotherParticleEvent[[2]][[2;;4]]];
eps=-norm*{0.085,0.085,0.085};

If[AxisBool===1,axis=Show[Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,0,norm*0.4}}]}],Graphics3D[Text["z",{norm*0.055,norm*0.055,norm*0.425}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,norm*0.4,0}}]}],Graphics3D[Text["y",{norm*0.055,norm*0.425,norm*0.055}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{norm*0.4,0,0}}]}],Graphics3D[Text["x",{norm*0.425,norm*0.055,norm*0.055}]]](*axis*);,axis=Graphics3D[]];


vertex2Visualize=Show[Graphics3D[{PointSize[0.012],Black,Point[VertexT]}],Graphics3D[Text[Style["\!\(\*SubscriptBox[\(v\), \(2\)]\)",Black,FontSize->14],VertexT+eps]]](*vertex decay 3-body*);

If[HideMotherMomenumFlag===1,motherVisualize=Show[Graphics3D[{Purple,Dashed,Thickness[0.006],Line[{MotherParticleEvent[[3]],VertexT}]}]](*mother particle*);,motherVisualize=Show[Graphics3D[{Purple,Arrowheads[0.03],Arrow[{MotherParticleEvent[[3]],MotherParticleEvent[[2]][[2;;4]]+MotherParticleEvent[[3]]}]}],Graphics3D[Text[Style[MotherParticleEvent[[1]],Purple,FontSize->14],MotherParticleEvent[[2]][[2;;4]]+MotherParticleEvent[[3]]+eps]],Graphics3D[{Purple,Dashed,Thickness[0.006],Line[{MotherParticleEvent[[3]],VertexT}]}]](*mother particle*);];

daugther1Visualize=Show[Graphics3D[{Purple,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[1]][[1]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[1]][[1]][[1]],Purple,FontSize->14],DaughtersEvent[[1]][[1]][[2]][[2;;4]]+VertexT+eps]]](*particle 1*);

daugther2Visualize=Show[Graphics3D[{Red,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[1]][[2]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[1]][[2]][[1]],Red,FontSize->14],DaughtersEvent[[1]][[2]][[2]][[2;;4]]+VertexT+eps]]](*particle 2*);

daugther3Visualize=Show[Graphics3D[{Green,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[1]][[3]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[1]][[3]][[1]],Green,FontSize->14],DaughtersEvent[[1]][[3]][[2]][[2;;4]]+VertexT+eps]]](*particle 3*);

Show[motherVisualize,vertex2Visualize,daugther1Visualize,daugther2Visualize,daugther3Visualize,axis,ImageSize->Large]
 ]


(* ::Subsubsection::Closed:: *)
(*Visualization ThreeBodySM*)


ThreeBodyEventsDisplaySM[MotherParticleDecayVertex_,DaughtersEvent_,AxisBool_]:=Module[{VertexT,norm,vertex2Visualize,eps,axis,motherVisualize,daugther1Visualize,daugther2Visualize,daugther3Visualize},

VertexT=MotherParticleDecayVertex;
norm=0.3*Norm[DaughtersEvent[[1]][[2]][[2;;4]]];
eps=-norm*{0.085,0.085,0.085};

If[AxisBool===1,axis=Show[Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,0,norm*0.4}}]}],Graphics3D[Text["z",{norm*0.055,norm*0.055,norm*0.425}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,norm*0.4,0}}]}],Graphics3D[Text["y",{norm*0.055,norm*0.425,norm*0.055}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{norm*0.4,0,0}}]}],Graphics3D[Text["x",{norm*0.425,norm*0.055,norm*0.055}]]](*axis*);,axis=Graphics3D[]];


vertex2Visualize=Show[Graphics3D[{PointSize[0.012],Black,Point[VertexT]}],Graphics3D[Text[Style["\!\(\*SubscriptBox[\(v\), \(2\)]\)",Black,FontSize->14],VertexT+eps]]](*vertex decay 3-body*);

daugther1Visualize=Show[Graphics3D[{Red,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[1]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[1]][[1]],Red,FontSize->14],DaughtersEvent[[1]][[2]][[2;;4]]+VertexT+eps]]](*particle 1*);

daugther2Visualize=Show[Graphics3D[{Green,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[2]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[2]][[1]],Green,FontSize->14],DaughtersEvent[[2]][[2]][[2;;4]]+VertexT+eps]]](*particle 2*);

Show[vertex2Visualize,daugther1Visualize,daugther2Visualize,axis,ImageSize->Large]
 ]


(* ::Subsubsection::Closed:: *)
(*Visualization TwoBody*)


TwoBodyEventsDisplay[MotherParticleEvent_,DaughtersEvent_,AxisBool_,HideMotherMomenumFlag_]:=Module[{VertexT,vertex1Visualize,eps,norm,axis,motherVisualize,daugther1Visualize,daugther2Visualize},

VertexT=DaughtersEvent[[2]];
norm=0.2*Norm[MotherParticleEvent[[2]][[2;;4]]];
eps=-norm*{0.085,0.085,0.085};

If[AxisBool===1,axis=Show[Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,0,norm*0.4}}]}],Graphics3D[Text["z",{norm*0.055,norm*0.055,norm*0.425}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,norm*0.4,0}}]}],Graphics3D[Text["y",{norm*0.055,norm*0.425,norm*0.055}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{norm*0.4,0,0}}]}],Graphics3D[Text["x",{norm*0.425,norm*0.055,norm*0.055}]]](*axis*);,axis=Graphics3D[]];


vertex1Visualize=Show[Graphics3D[{PointSize[0.012],Black,Point[VertexT]}],Graphics3D[Text[Style["\!\(\*SubscriptBox[\(v\), \(1\)]\)",Black,FontSize->14],VertexT+eps]]](*vertex decay 2-body*);


If[HideMotherMomenumFlag===1,motherVisualize=Show[Graphics3D[{Black,Dashed,Thickness[0.006],Line[{MotherParticleEvent[[3]],VertexT}]}]](*mother particle*);,motherVisualize=Show[Graphics3D[{Black,Arrowheads[0.03],Arrow[{MotherParticleEvent[[3]],MotherParticleEvent[[2]][[2;;4]]+MotherParticleEvent[[3]]}]}],Graphics3D[Text[Style[MotherParticleEvent[[1]],Black,FontSize->14],MotherParticleEvent[[2]][[2;;4]]+MotherParticleEvent[[3]]+eps]],Graphics3D[{Black,Dashed,Thickness[0.006],Line[{MotherParticleEvent[[3]],VertexT}]}]](*mother particle*);];


daugther1Visualize=Show[Graphics3D[{Purple,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[1]][[1]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[1]][[1]][[1]],Purple,FontSize->14],DaughtersEvent[[1]][[1]][[2]][[2;;4]]+VertexT+eps]]](*particle 1*);

daugther2Visualize=Show[Graphics3D[{Purple,Arrowheads[0.03],Arrow[{VertexT,DaughtersEvent[[1]][[2]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[DaughtersEvent[[1]][[2]][[1]],Purple,FontSize->14],DaughtersEvent[[1]][[2]][[2]][[2;;4]]+VertexT+eps]]](*particle 2*);

Show[motherVisualize,vertex1Visualize,daugther1Visualize,daugther2Visualize,axis,ImageSize->Large]
 ]


(* ::Subsubsection::Closed:: *)
(*Visualization FullEvent - assumes at least one ch-*)


FullEventsDisplay["e- e+ -> \[Gamma] Z'",FullEvent_,AxisBool_,HidePhotonFlag_]:=Module[{VertexT,vertex0Visualize,eps,norm,axis,photonVisualize,ZpDecayVisualize,ChiPlusDecacyVisualize},

VertexT={0,0,0};
norm=0.4*Norm[FullEvent[[2]][[1]][[2]][[2;;4]]];
eps=-norm*{0.05,0.05,0.05};

If[AxisBool===1,axis=Show[Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,0,norm*0.4}}]}],Graphics3D[Text["z",{norm*0.05,norm*0.05,norm*0.4}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,norm*0.4,0}}]}],Graphics3D[Text["y",{norm*0.05,norm*0.4,norm*0.05}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{norm*0.4,0,0}}]}],Graphics3D[Text["x",{norm*0.4,norm*0.05,norm*0.05}]]](*axis*);,axis=Graphics3D[]];


vertex0Visualize=Show[Graphics3D[{PointSize[0.012],Black,Point[VertexT]}],Graphics3D[Text[Style["\!\(\*SubscriptBox[\(v\), \(0\)]\)",Black,FontSize->14],VertexT+eps]]](*vertex of collision*);

If[HidePhotonFlag===1,photonVisualize=Graphics3D[];,
photonVisualize=Show[Graphics3D[{Blue,Arrowheads[0.03],Arrow[{VertexT,FullEvent[[1]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[FullEvent[[1]][[1]],Blue,FontSize->14],FullEvent[[1]][[2]][[2;;4]]+VertexT+eps]]](*photon*);];

ZpDecayVisualize= TwoBodyEventsDisplay[{FullEvent[[2]][[1]][[1]],FullEvent[[2]][[1]][[2]],VertexT},{{FullEvent[[2]][[2]][[1]][[1]],FullEvent[[2]][[2]][[1]][[2]][[1]]},FullEvent[[2]][[2]][[2]]},0,0];

ChiPlusDecacyVisualize= ThreeBodyEventsDisplay[{FullEvent[[2]][[2]][[1]][[2]][[1]][[1]],FullEvent[[2]][[2]][[1]][[2]][[1]][[2]],FullEvent[[2]][[2]][[2]]},FullEvent[[2]][[2]][[1]][[2]][[2]],0,1];

Show[vertex0Visualize,photonVisualize,ZpDecayVisualize,ChiPlusDecacyVisualize,axis,ImageSize->Large,ViewVertical->{1,0,0},ViewPoint->{1.5,-2.5,-1.25}]
 ]


(* ::Subsubsection::Closed:: *)
(*Visualization SMEvent - assumes at least one ch-*)


SMEventsDisplay["e- e+ -> \[Gamma] Z'",SMEvent_,AxisBool_]:=Module[{VertexT,vertex0Visualize,eps,norm,axis,photonVisualize,ZpDecayVisualize,ChiPlusDecacyVisualizeSM},

VertexT={0,0,0};
norm=0.4*Norm[SMEvent[[1]][[2]][[2;;4]]];
eps=-norm*{0.05,0.05,0.05};

If[AxisBool===1,axis=Show[Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,0,norm*0.4}}]}],Graphics3D[Text["z",{norm*0.05,norm*0.05,norm*0.4}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{0,norm*0.4,0}}]}],Graphics3D[Text["y",{norm*0.05,norm*0.4,norm*0.05}]],Graphics3D[{Opacity[0.2,Black],Arrowheads[0.02],Arrow[{{0,0,0},{norm*0.4,0,0}}]}],Graphics3D[Text["x",{norm*0.4,norm*0.05,norm*0.05}]]](*axis*);,axis=Graphics3D[]];


vertex0Visualize=Show[Graphics3D[{PointSize[0.012],Black,Point[VertexT]}],Graphics3D[Text[Style["\!\(\*SubscriptBox[\(v\), \(0\)]\)",Black,FontSize->14],VertexT+eps]]](*vertex of collision*);

photonVisualize=Show[Graphics3D[{Blue,Arrowheads[0.03],Arrow[{VertexT,SMEvent[[1]][[2]][[2;;4]]+VertexT}]}],Graphics3D[Text[Style[SMEvent[[1]][[1]],Blue,FontSize->14],SMEvent[[1]][[2]][[2;;4]]+VertexT+eps]]](*photon*);

ChiPlusDecacyVisualizeSM= ThreeBodyEventsDisplaySM[SMEvent[[3]][[1]],{SMEvent[[3]][[2]],SMEvent[[3]][[3]]},0];


Show[vertex0Visualize,photonVisualize,ChiPlusDecacyVisualizeSM,axis,ImageSize->Large,ViewVertical->{1,0,0},ViewPoint->{1.5,-2.5,-1.25}]
 ]


(* ::Chapter::Closed:: *)
(*Modules for compilation*)


(* ::Section::Closed:: *)
(*Boosts*)


DoubleFourMomentaBoost=Null;
TripleFourMomentaBoost=Null;
CompileBoosts[]:=Module[
{duration,(*functions*)EBoostedByV,PxBoostedByV,PyBoostedByV,PzBoostedByV,(*variables*)Ener,px,py,pz,vx,vy,vz,(*variables for compilation*)temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,E1,p1x,p1y,p1z,E2,p2x,p2y,p2z,E3,p3x,p3y,p3z},

duration=Now;

(*functions*)
EBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=Ener/Sqrt[1-vx^2-vy^2-vz^2]-(px*vx)/Sqrt[1-vx^2-vy^2-vz^2]-(py*vy)/Sqrt[1-vx^2-vy^2-vz^2]-(pz*vz)/Sqrt[1-vx^2-vy^2-vz^2];
PxBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=-((Ener*vx)/Sqrt[1-vx^2-vy^2-vz^2])+(py*vx*vy*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+(pz*vx*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+px*(1+(vx^2*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2));
PyBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=-((Ener*vy)/Sqrt[1-vx^2-vy^2-vz^2])+(px*vx*vy*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+(pz*vy*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+py*(1+(vy^2*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2));
PzBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=-((Ener*vz)/Sqrt[1-vx^2-vy^2-vz^2])+(px*vx*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+(py*vy*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+pz*(1+(vz^2*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2));


(*variables for compilation*)
temp1=EBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];
temp2=PxBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];
temp3=PyBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];
temp4=PzBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];

temp5=EBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];
temp6=PxBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];
temp7=PyBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];
temp8=PzBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];

temp9=EBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];
temp10=PxBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];
temp11=PyBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];
temp12=PzBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];


(*compilation - functions already defined outside*)
DoubleFourMomentaBoost=
Compile[{
{E1, _Real},{p1x,_Real},{p1y,_Real},{p1z,_Real},{E2, _Real},{p2x,_Real},{p2y,_Real},{p2z,_Real},{vx, _Real},{vy,_Real},{vz,_Real}},
{{temp1,temp2,temp3,temp4},{temp5,temp6,temp7,temp8}},
CompilationTarget->"C",
Parallelization->True,
RuntimeAttributes->{Listable},RuntimeOptions->"Speed",
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];

TripleFourMomentaBoost=
Compile[{
{E1, _Real},{p1x,_Real},{p1y,_Real},{p1z,_Real},{E2, _Real},{p2x,_Real},{p2y,_Real},{p2z,_Real},{E3, _Real},{p3x,_Real},{p3y,_Real},{p3z,_Real},{vx, _Real},{vy,_Real},{vz,_Real}},
{{temp1,temp2,temp3,temp4},{temp5,temp6,temp7,temp8},{temp9,temp10,temp11,temp12}},
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
(*RuntimeAttributes\[Rule]{Listable},*)
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];


(*checking if functions were compiled - if not, just return a usual mathematica function*)
If[(StringContainsQ[CompilePrint[DoubleFourMomentaBoost],"MainEvaluate"])||(StringContainsQ[CompilePrint[TripleFourMomentaBoost],"MainEvaluate"]),Print["ERROR: boosts not compiled. Proceeding with uncompiled boosts."];DoubleFourMomentaBoost[E1list_,p1xlist_,p1ylist_,p1zlist_,E2list_,p2xlist_,p2ylist_,p2zlist_,vxlist_,vylist_,vzlist_]:=Module[{i},Table[{{temp1,temp2,temp3,temp4},{temp5,temp6,temp7,temp8}}/.{E1->E1list[[i]],p1x->p1xlist[[i]],p1y->p1ylist[[i]],p1z->p1zlist[[i]],E2->E2list[[i]],p2x->p2xlist[[i]],p2y->p2ylist[[i]],p2z->p2zlist[[i]],vx->vxlist[[i]],vy->vylist[[i]],vz->vzlist[[i]]},{i,1,E1list//Length}]];TripleFourMomentaBoost[E1list_,p1xlist_,p1ylist_,p1zlist_,E2list_,p2xlist_,p2ylist_,p2zlist_,E3list_,p3xlist_,p3ylist_,p3zlist_,vxlist_,vylist_,vzlist_]:=Module[{i},Table[{{temp1,temp2,temp3,temp4},{temp5,temp6,temp7,temp8},{temp9,temp10,temp11,temp12}}/.{E1->E1list[[i]],p1x->p1xlist[[i]],p1y->p1ylist[[i]],p1z->p1zlist[[i]],E2->E2list[[i]],p2x->p2xlist[[i]],p2y->p2ylist[[i]],p2z->p2zlist[[i]],E3->E3list[[i]],p3x->p3xlist[[i]],p3y->p3ylist[[i]],p3z->p3zlist[[i]],vx->vxlist[[i]],vy->vylist[[i]],vz->vzlist[[i]]},{i,1,E1list//Length}]];,Print["\[Checkmark] Boosts compiled."];];

duration=(Now-duration)//First
]


CompileBoosts[0]:=Module[
{duration,(*functions*)EBoostedByV,PxBoostedByV,PyBoostedByV,PzBoostedByV,(*variables*)Ener,px,py,pz,vx,vy,vz,temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9,temp10,temp11,temp12,E1,p1x,p1y,p1z,E2,p2x,p2y,p2z,E3,p3x,p3y,p3z},

duration=Now;

(*functions*)
EBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=Ener/Sqrt[1-vx^2-vy^2-vz^2]-(px*vx)/Sqrt[1-vx^2-vy^2-vz^2]-(py*vy)/Sqrt[1-vx^2-vy^2-vz^2]-(pz*vz)/Sqrt[1-vx^2-vy^2-vz^2];
PxBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=-((Ener*vx)/Sqrt[1-vx^2-vy^2-vz^2])+(py*vx*vy*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+(pz*vx*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+px*(1+(vx^2*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2));
PyBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=-((Ener*vy)/Sqrt[1-vx^2-vy^2-vz^2])+(px*vx*vy*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+(pz*vy*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+py*(1+(vy^2*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2));
PzBoostedByV[Ener_,px_,py_,pz_,vx_,vy_,vz_]=-((Ener*vz)/Sqrt[1-vx^2-vy^2-vz^2])+(px*vx*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+(py*vy*vz*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2)+pz*(1+(vz^2*(-1+1/Sqrt[1-vx^2-vy^2-vz^2]))/(vx^2+vy^2+vz^2));


(*variables for compilation*)
temp1=EBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];
temp2=PxBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];
temp3=PyBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];
temp4=PzBoostedByV[E1,p1x,p1y,p1z,vx,vy,vz];

temp5=EBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];
temp6=PxBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];
temp7=PyBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];
temp8=PzBoostedByV[E2,p2x,p2y,p2z,vx,vy,vz];

temp9=EBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];
temp10=PxBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];
temp11=PyBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];
temp12=PzBoostedByV[E3,p3x,p3y,p3z,vx,vy,vz];


(*checking if functions were compiled - if not, just return a usual mathematica function*)
Print["Proceeding with uncompiled boosts"];DoubleFourMomentaBoost[E1list_,p1xlist_,p1ylist_,p1zlist_,E2list_,p2xlist_,p2ylist_,p2zlist_,vxlist_,vylist_,vzlist_]:=Module[{i},Table[{{temp1,temp2,temp3,temp4},{temp5,temp6,temp7,temp8}}/.{E1->E1list[[i]],p1x->p1xlist[[i]],p1y->p1ylist[[i]],p1z->p1zlist[[i]],E2->E2list[[i]],p2x->p2xlist[[i]],p2y->p2ylist[[i]],p2z->p2zlist[[i]],vx->vxlist[[i]],vy->vylist[[i]],vz->vzlist[[i]]},{i,1,E1list//Length}]];TripleFourMomentaBoost[E1list_,p1xlist_,p1ylist_,p1zlist_,E2list_,p2xlist_,p2ylist_,p2zlist_,E3list_,p3xlist_,p3ylist_,p3zlist_,vxlist_,vylist_,vzlist_]:=Module[{i},Table[{{temp1,temp2,temp3,temp4},{temp5,temp6,temp7,temp8},{temp9,temp10,temp11,temp12}}/.{E1->E1list[[i]],p1x->p1xlist[[i]],p1y->p1ylist[[i]],p1z->p1zlist[[i]],E2->E2list[[i]],p2x->p2xlist[[i]],p2y->p2ylist[[i]],p2z->p2zlist[[i]],E3->E3list[[i]],p3x->p3xlist[[i]],p3y->p3ylist[[i]],p3z->p3zlist[[i]],vx->vxlist[[i]],vy->vylist[[i]],vz->vzlist[[i]]},{i,1,E1list//Length}]];


duration=(Now-duration)//First
]


(* ::Section::Closed:: *)
(*Boundaries for the Dalitz plane*)


m23minCompiled=Null;
m23maxCompiled=Null;
CompileDalitzBoundaries[]:=Module[
{duration,(*functions*)m23min,m23max,(*variables*)m12,M,m1,m2,m3,(*variables for compilation*)temp1,temp2},

duration=Now;

(*functions*)
m23min[m12_,M_,m1_,m2_,m3_]=1/(2m12) (M^2 (-m1^2+m12+m2^2)-(m12-m2^2) (m12-m3^2)+m1^2 (m12+m3^2)-Sqrt[(m1^2-m12)^2-2 (m1^2+m12) m2^2+m2^4] Sqrt[(M^2-m12)^2-2 (M^2+m12) m3^2+m3^4]);
m23max[m12_,M_,m1_,m2_,m3_]=1/(2m12) (M^2 (-m1^2+m12+m2^2)-(m12-m2^2) (m12-m3^2)+m1^2 (m12+m3^2)+Sqrt[(m1^2-m12)^2-2 (m1^2+m12) m2^2+m2^4] Sqrt[(M^2-m12)^2-2 (M^2+m12) m3^2+m3^4]);


(*variables for compilation*)
temp1=m23min[m12,M,m1,m2,m3];
temp2=m23max[m12,M,m1,m2,m3];


(*compilation - functions already defined outside*)
m23minCompiled=
Compile[{
{M, _Real},{m1,_Real},{m2,_Real},{m3,_Real},{m12,_Real}},
temp1,
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
(*RuntimeAttributes\[Rule]{Listable},*)
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];

m23maxCompiled=
Compile[{
{M, _Real},{m1,_Real},{m2,_Real},{m3,_Real},{m12,_Real}},
temp2,
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
(*RuntimeAttributes\[Rule]{Listable},*)
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];


(*checking if functions were compiled - if not, just return a usual mathematica function*)
If[(StringContainsQ[CompilePrint[m23minCompiled],"MainEvaluate"])||(StringContainsQ[CompilePrint[m23maxCompiled],"MainEvaluate"]),Print["ERROR: boundaries for the Dalitz plane not compiled. Proceeding with uncompiled boundaries."];m23minCompiled=m23min;m23maxCompiled=m23max;,Print["\[Checkmark] Boundaries for the Dalitz plane compiled."];];

duration=(Now-duration)//First
]


CompileDalitzBoundaries[0]:=Module[
{duration,(*functions*)m23min,m23max,(*variables*)m12,M,m1,m2,m3,(*variables for compilation*)temp1,temp2},

duration=Now;

(*functions*)
m23min[m12_,M_,m1_,m2_,m3_]=1/(2m12) (M^2 (-m1^2+m12+m2^2)-(m12-m2^2) (m12-m3^2)+m1^2 (m12+m3^2)-Sqrt[(m1^2-m12)^2-2 (m1^2+m12) m2^2+m2^4] Sqrt[(M^2-m12)^2-2 (M^2+m12) m3^2+m3^4]);
m23max[m12_,M_,m1_,m2_,m3_]=1/(2m12) (M^2 (-m1^2+m12+m2^2)-(m12-m2^2) (m12-m3^2)+m1^2 (m12+m3^2)+Sqrt[(m1^2-m12)^2-2 (m1^2+m12) m2^2+m2^4] Sqrt[(M^2-m12)^2-2 (M^2+m12) m3^2+m3^4]);


Print["Proceeding with uncompiled boundaries."];m23minCompiled=m23min;m23maxCompiled=m23max;

duration=(Now-duration)//First
]


(* ::Section::Closed:: *)
(*Electron-positron collision*)


ElectronPositronCollisionDistributionNoEpsCompiled=Null;
AssigningMomentaElectronPositronCollisions=Null;
ReturningParticlesElectronPositronCollisions=Null;
CompileElectronPositronCollisions[]:=Module[
{duration,(*functions*)t0,ElectronPositronCollisionMandelstam,ElectronPositronCollisionDistribution,(*variables*)mZp,s,t,thetaCM,(*variables for compilation*)temp1,(*fixed values*)vME},

duration=Now;

(*definitions - fixed values*)
vME=ElectronMass;


(*functions*)
t0[mZp_,s_]=1/2 (Sqrt[(mZp^2-s)^2/s] Sqrt[s-4 vME^2]+mZp^2-s+2 vME^2);
ElectronPositronCollisionMandelstam[s_,t_,mZp_]=-(1/((t-vME^2)^2 (mZp^2-s-t+vME^2)^2))2   (-(mZp^6 (t-3 vME^2))+mZp^4 (s (t-5 vME^2)+3 t^2-6 t vME^2+7 vME^4)-mZp^2 (s^2 (t-3 vME^2)+4 s (t^2-t vME^2+2 vME^4)+4 (t-vME^2)^3)+s^3 (t-vME^2)+s^2 (3 t^2-2 t vME^2+3 vME^4)+4 s t (t-vME^2)^2+2 (t-vME^2)^4);
ElectronPositronCollisionDistribution[s_,thetaCM_,mZp_]=ElectronPositronCollisionMandelstam[s,t0[mZp,s]-4Sqrt[1/4 (s-4 vME^2)] Sqrt[(mZp^2-s)^2/(4 s)] Sin[thetaCM/2]^2,mZp](*no need for epsilon here since we just use this as a distribution, so it is set to 1*);


(*variables for compilation*)
temp1=ElectronPositronCollisionDistribution[s,thetaCM,mZp];


(*compilation - functions already defined outside*)
ElectronPositronCollisionDistributionNoEpsCompiled=
Compile[{
{s, _Real},{thetaCM,_Real},{mZp,_Real}},
temp1,
CompilationTarget->"C",
Parallelization->True,
RuntimeAttributes->{Listable},RuntimeOptions->"Speed",
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];

AssigningMomentaElectronPositronCollisions=Compile[{
{E1, _Real},{p1, _Real},{E2, _Real},{vx, _Real},{vy, _Real},{vz, _Real},{phiCMs, _Real},{thetaCMs,_Real},{nEvents,_Integer}},
{E1,p1*Cos[phiCMs]*Sin[thetaCMs],p1*Sin[phiCMs]*Sin[thetaCMs],p1*Cos[thetaCMs],E2,-p1*Cos[phiCMs]*Sin[thetaCMs],-p1*Sin[phiCMs]*Sin[thetaCMs],-p1*Cos[thetaCMs],-vx,-vy,-vz},
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
RuntimeAttributes->{Listable},
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}] ;

ReturningParticlesElectronPositronCollisions=Compile[{
{FinalStateParticles1,_Integer},{FinalStateParticles2,_Integer},{FourMomentaLAB,_Real,2},{nEvents,_Integer}},
{{{FinalStateParticles1,FinalStateParticles1,FinalStateParticles1,FinalStateParticles1},{Compile`GetElement[FourMomentaLAB,1,1],Compile`GetElement[FourMomentaLAB,1,2],Compile`GetElement[FourMomentaLAB,1,3],Compile`GetElement[FourMomentaLAB,1,4]}},{{FinalStateParticles2,FinalStateParticles2,FinalStateParticles2,FinalStateParticles2},{Compile`GetElement[FourMomentaLAB,2,1],Compile`GetElement[FourMomentaLAB,2,2],Compile`GetElement[FourMomentaLAB,2,3],Compile`GetElement[FourMomentaLAB,2,4]}}},
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
RuntimeAttributes->{Listable},
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}] ;


(*checking if functions were compiled - if not, just return a usual mathematica function*)
If[(StringContainsQ[CompilePrint[ElectronPositronCollisionDistributionNoEpsCompiled],"MainEvaluate"])||(StringContainsQ[CompilePrint[ElectronPositronCollisionDistributionNoEpsCompiled],"MainEvaluate"])||(StringContainsQ[CompilePrint[ElectronPositronCollisionDistributionNoEpsCompiled],"MainEvaluate"]),Print["ERROR: electron-positron collision distribution not compiled. Proceeding with uncompiled distribution."];ElectronPositronCollisionDistributionNoEpsCompiled[s_,thetaCM_,mZp_]:=Module[{i},Table[ElectronPositronCollisionMandelstam[s,thetaCM[[i]],mZp],{i,1,thetaCM//Length}]];
AssigningMomentaElectronPositronCollisions[E1_,p1_,E2_,vx_,vy_,vz_,\[Phi]CMs_,\[Theta]CMs_,nEvents_]:=Module[{px,py,pz,i},Table[{E1,px=p1*Cos[\[Phi]CMs[[i]]]*Sin[\[Theta]CMs[[i]]],py=p1*Sin[\[Phi]CMs[[i]]]*Sin[\[Theta]CMs[[i]]],pz=p1*Cos[\[Theta]CMs[[i]]],E2,-px,-py,-pz,-vx,-vy,-vz},{i,1,nEvents}]];
ReturningParticlesElectronPositronCollisions[FinalStateParticles1_,FinalStateParticles2_,FourMomentaLAB_,nEvents_]:=Module[{i},Table[{{FinalStateParticles1,FourMomentaLAB[[i]][[1]]},{FinalStateParticles2,FourMomentaLAB[[i]][[2]]}},{i,1,nEvents}]];,Print["\[Checkmark] Electron-positron collision distribution compiled."];];

duration=(Now-duration)//First
]


CompileElectronPositronCollisions[0]:=Module[
{duration,(*functions*)t0,ElectronPositronCollisionMandelstam,ElectronPositronCollisionDistribution,(*variables*)mZp,s,t,thetaCM,(*variables for compilation*)temp1,(*fixex values*)vME},

duration=Now;

(*fixex values*)
vME=0.5109989461*10^-3(*GeV*);


(*functions*)
t0[mZp_,s_]=1/2 (Sqrt[(mZp^2-s)^2/s] Sqrt[s-4 vME^2]+mZp^2-s+2 vME^2);
ElectronPositronCollisionMandelstam[s_,t_,mZp_]=-(1/((t-vME^2)^2 (mZp^2-s-t+vME^2)^2))2   (-(mZp^6 (t-3 vME^2))+mZp^4 (s (t-5 vME^2)+3 t^2-6 t vME^2+7 vME^4)-mZp^2 (s^2 (t-3 vME^2)+4 s (t^2-t vME^2+2 vME^4)+4 (t-vME^2)^3)+s^3 (t-vME^2)+s^2 (3 t^2-2 t vME^2+3 vME^4)+4 s t (t-vME^2)^2+2 (t-vME^2)^4);
ElectronPositronCollisionDistribution[s_,thetaCM_,mZp_]=ElectronPositronCollisionMandelstam[s,t0[mZp,s]-4Sqrt[1/4 (s-4 vME^2)] Sqrt[(mZp^2-s)^2/(4 s)] Sin[thetaCM/2]^2,mZp](*no need for epsilon here since we just use this as a distribution, so it is set to 1*);



Print["Proceeding with uncompiled distribution."];ElectronPositronCollisionDistributionNoEpsCompiled[s_,thetaCM_,mZp_]:=Module[{i},Table[ElectronPositronCollisionMandelstam[s,thetaCM[[i]],mZp],{i,1,thetaCM//Length}]];
AssigningMomentaElectronPositronCollisions[E1_,p1_,E2_,vx_,vy_,vz_,\[Phi]CMs_,\[Theta]CMs_,nEvents_]:=Module[{px,py,pz,i},Table[{E1,px=p1*Cos[\[Phi]CMs[[i]]]*Sin[\[Theta]CMs[[i]]],py=p1*Sin[\[Phi]CMs[[i]]]*Sin[\[Theta]CMs[[i]]],pz=p1*Cos[\[Theta]CMs[[i]]],E2,-px,-py,-pz,-vx,-vy,-vz},{i,1,nEvents}]];
ReturningParticlesElectronPositronCollisions[FinalStateParticles1_,FinalStateParticles2_,FourMomentaLAB_,nEvents_]:=Module[{i},Table[{{FinalStateParticles1,FourMomentaLAB[[i]][[1]]},{FinalStateParticles2,FourMomentaLAB[[i]][[2]]}},{i,1,nEvents}]];

duration=(Now-duration)//First
]


(* ::Section::Closed:: *)
(*Chi+ decay*)


ChiDistributionElectrons3bodyCompiled=Null;
ChiDistributionMuons3bodyCompiled=Null;
ChiDistributionPions3bodyCompiled=Null;
CompileChiDecays[]:=Module[
{duration,(*functions*)ChiDistributionLeptons3body,ChiDistributionPions3body,(*variables*)M,m1,m2,m3,m12,m23,mZp,(*variables for compilation*)temp1,temp2,temp3},

duration=Now;


(*functions*)
ChiDistributionLeptons3body[M_,m1_,m2_,m3_,mZp_,m12_,m23_]=1/(mZp^4 (m23-mZp^2)^2) 1.7231563788153554` (1/(mZp^2-8315.178393760001`) 0.12132142105910404` ((2 m1^2+m2^2+m3^2-2 m12-m23-2 m2 m3) M^2-2 m1 (m2^2-4 m3 m2+m3^2-m23) M+2 m12^2-2 m12 m2^2+m23^2+2 m2^2 m3^2-2 m12 m3^2-m23 m3^2+m1^2 (m2^2-2 m3 m2+m3^2-2 m12-m23)-m2^2 m23+2 m12 m23+2 m2 m23 m3) mZp^6-2.866054752661671` ((2 m1^2+m2^2+m3^2-2 m12-m23-2 m2 m3) M^2-2 m1 (m2^2-4 m3 m2+m3^2-m23) M+2 m12^2-2 m12 m2^2+m23^2+2 m2^2 m3^2-2 m12 m3^2-m23 m3^2+m1^2 (m2^2-2 m3 m2+m3^2-2 m12-m23)-m2^2 m23+2 m12 m23+2 m2 m23 m3) mZp^4-4 (M-m1) (-((0.06066071052955202` (m2-m3) (m2 M^3+m1 m2 M^2+(m2^3+m3 m2^2-m23 m2+m1^2 m3-m12 (m2+m3)) M+m1 (m3 (m1^2-m23+m3 (m2+m3))-m12 (m2+m3))) mZp^2)/(mZp^2-8315.178393760001`))+1.4330273763308354` (m2-m3) (m2 M^3+m1 m2 M^2+(m2^3+m3 m2^2-m23 m2+m1^2 m3-m12 (m2+m3)) M+m1 (m3 (m1^2-m23+m3 (m2+m3))-m12 (m2+m3)))+1/(mZp^2-8315.178393760001`)^2 0.30511302550508157` (0.5021039710465376` M m2^4-0.5021039710465376` M (m3^2+m23) m2^2+M^3 (0.5021039710465376` m2+0.49789602895346247` m3) m2+M^2 m1 (0.5021039710465376` m2+0.49789602895346247` m3) m2-0.49789602895346247` M (m1^2+m23) m3 m2-0.5021039710465376` M m1^2 m3^2-0.5021039710465376` M m12 (m2^2-1.` m3^2)-0.5021039710465376` m1 m12 (m2^2-1.` m3^2)-0.5021039710465376` m1 m3 (1.` m3^3-1.` m2^2 m3-1.` m23 m3-0.9916193809734177` m2 m23+m1^2 (0.9916193809734177` m2+1.` m3))) mZp^4) mZp^2+(M-m1)^2 (1/(mZp^2-8315.178393760001`)^2 0.30511302550508157` (-0.5021039710465376` m2^4+0.5021039710465376` (2 m3^2+m23) m2^2+0.995792057906925` m23 m3 m2-0.5021039710465376` m3^2 (m3^2-1.` m23)) mZp^4+(0.06066071052955202` (m2^2+2 m3 m2+m3^2-m23) (m2-m3)^2 mZp^2)/(mZp^2-8315.178393760001`)- 1.4330273763308354` (m2-m3)^2 (m2^2+2 m3 m2+m3^2-m23)) (M^2+2 m1 M+m1^2-m23)-1/(mZp^2-8315.178393760001`)^2 0.6127938468964998` ((1.` m1^2+0.5` m2^2+0.5` m3^2-1.` m12-0.4999999999999999` m23+0.9916193809734177` m2 m3) M^2+m1 (-1.` m2^2-3.9664775238936714` m3 m2-1.` m3^2+1.` m23) M+1.` m12^2-0.9999999999999998` m12 m2^2+0.5` m23^2+1.` m2^2 m3^2-0.9999999999999998` m12 m3^2-0.4999999999999999` m23 m3^2+m1^2 (0.5` m2^2+0.9916193809734178` m3 m2+0.5` m3^2-1.` m12-0.5` m23)-0.4999999999999999` m2^2 m23+1.` m12 m23-0.9916193809734177` m2 m23 m3) mZp^8) ;
ChiDistributionPions3body[M_,m1_,m2_,m3_,mZp_,m12_,m23_]=1/(mZp^4 (m23-mZp^2)^2) 0.6173325661353781`  (1.1265731376593378` Abs[((0.6547926192842074` +0.00001179288147985907` I) (1+(0.001004506501192573` -0.000011050017200740626` I)/((0.6125410225` -0.006644698499999999` I)-m23)))/((0.-0.23479121194530173` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+1.8783296955624138` (0.04271767567565097` (0.6010280675999999` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.4990469642661888`))+0.6010280675999999`)+(0.020154984130451314` +0.0169763089414969` I)/((0.-0.04949270045567086` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+0.3959416036453669` (0.040860712191285925` (5.080516` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.8774686004328772`))+5.080516`)+(0.044963937180827813` +0.24598370047761595` I)/((0.-0.17569710967673408` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+1.4055768774138726` (0.041181101512326464` (3.463321` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.8134044037295003`))+3.463321`)-(0.3183349993211615` +0.22649769503854025` I)/((0.-0.3016808572181247` I) m23 (1-0.07791982070308841`/m23)^(3/2)-m23+2.4134468577449977` (0.041623970301305685` (2.2290490000000003` -m23)+1/4 (1-0.07791982070308841`/m23) m23 (1/\[Pi] Sqrt[1-0.07791982070308841`/m23] Log[3.5824160974864263` (Sqrt[1-0.07791982070308841`/m23] Sqrt[m23]+Sqrt[m23])]-0.7383136955592318`))+2.2290490000000003`)])^2 (2 mZp^4 (M^4+M^2 (2 m1^2-4 m12-m23)+2 M m1 (2 m2^2-m23+2 m3^2)+m1^4-m1^2 (4 m12+m23)+4 (m12^2-m12 (m2^2-m23+m3^2)+m2^2 m3^2))+2 (M-m1) (m2^2-m3^2) ((M-m1) (m2^2-m3^2) (M^2+2 M m1+m1^2-m23)+2 mZp^2 (M^3+M^2 m1+M (m1^2-2 m12+2 m2^2-m23)+m1 (m1^2-2 m12-m23+2 m3^2))));


(*variables for compilation*)
temp1=ChiDistributionLeptons3body[M,m1,ElectronMass,ElectronMass,mZp,m12,m23](*only masses are important for the shape*);
temp2=ChiDistributionLeptons3body[M,m1,MuonMass,MuonMass,mZp,m12,m23](*only masses are important for the shape*);
temp3=ChiDistributionPions3body[M,m1,PionMass,PionMass,mZp,m12,m23](*only masses are important for the shape*);


(*compilation - functions already defined outside*)
ChiDistributionElectrons3bodyCompiled=
Compile[{
{M,_Real},{m1, _Real},{mZp, _Real},{m12, _Real},{m23, _Real}},
temp1,
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
(*RuntimeAttributes\[Rule]{Listable},*)
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];

ChiDistributionMuons3bodyCompiled=
Compile[{
{M,_Real},{m1, _Real},{mZp, _Real},{m12, _Real},{m23, _Real}},
temp2,
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
(*RuntimeAttributes\[Rule]{Listable},*)
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];

ChiDistributionPions3bodyCompiled=
Compile[{
{M,_Real},{m1, _Real},{mZp, _Real},{m12, _Real},{m23, _Real}},
temp3,
CompilationTarget->"C",
Parallelization->True,RuntimeOptions->"Speed",
(*RuntimeAttributes\[Rule]{Listable},*)
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];


(*checking if functions were compiled - if not, just return a usual mathematica function*)
If[(StringContainsQ[CompilePrint[ChiDistributionElectrons3bodyCompiled],"MainEvaluate"])||(StringContainsQ[CompilePrint[ChiDistributionMuons3bodyCompiled],"MainEvaluate"])||(StringContainsQ[CompilePrint[ChiDistributionPions3bodyCompiled],"MainEvaluate"]),Print["ERROR: chi+ decay distribution not compiled. Proceeding with uncompiled distribution."];ChiDistributionElectrons3bodyCompiled[M_,m1_,mZp_,m12_,m23_]=ChiDistributionLeptons3body[M,m1,ElectronMass,ElectronMass,mZp,m12,m23];
ChiDistributionMuons3bodyCompiled[M_,m1_,mZp_,m12_,m23_]=ChiDistributionLeptons3body[M,m1,MuonMass,MuonMass,mZp,m12,m23];ChiDistributionPions3bodyCompiled[M_,m1_,mZp_,m12_,m23_]=ChiDistributionPions3body[M,m1,PionMass,PionMass,mZp,m12,m23];,Print["\[Checkmark] Electron chi+ decay distribution compiled."];];

duration=(Now-duration)//First
]


CompileChiDecays[0]:=Module[
{duration,(*functions*)ChiDistributionLeptons3body,(*variables*)M,m1,m2,m3,m12,m23,mZp,(*variables for compilation*)temp1,temp2},

duration=Now;


(*functions*)
ChiDistributionLeptons3body[M_,m1_,m2_,m3_,mZp_,m12_,m23_]=1/(mZp^4 (m23-mZp^2)^2) 1.7231563788153554` (1/(mZp^2-8315.178393760001`) 0.12132142105910404` ((2 m1^2+m2^2+m3^2-2 m12-m23-2 m2 m3) M^2-2 m1 (m2^2-4 m3 m2+m3^2-m23) M+2 m12^2-2 m12 m2^2+m23^2+2 m2^2 m3^2-2 m12 m3^2-m23 m3^2+m1^2 (m2^2-2 m3 m2+m3^2-2 m12-m23)-m2^2 m23+2 m12 m23+2 m2 m23 m3) mZp^6-2.866054752661671` ((2 m1^2+m2^2+m3^2-2 m12-m23-2 m2 m3) M^2-2 m1 (m2^2-4 m3 m2+m3^2-m23) M+2 m12^2-2 m12 m2^2+m23^2+2 m2^2 m3^2-2 m12 m3^2-m23 m3^2+m1^2 (m2^2-2 m3 m2+m3^2-2 m12-m23)-m2^2 m23+2 m12 m23+2 m2 m23 m3) mZp^4-4 (M-m1) (-((0.06066071052955202` (m2-m3) (m2 M^3+m1 m2 M^2+(m2^3+m3 m2^2-m23 m2+m1^2 m3-m12 (m2+m3)) M+m1 (m3 (m1^2-m23+m3 (m2+m3))-m12 (m2+m3))) mZp^2)/(mZp^2-8315.178393760001`))+1.4330273763308354` (m2-m3) (m2 M^3+m1 m2 M^2+(m2^3+m3 m2^2-m23 m2+m1^2 m3-m12 (m2+m3)) M+m1 (m3 (m1^2-m23+m3 (m2+m3))-m12 (m2+m3)))+1/(mZp^2-8315.178393760001`)^2 0.30511302550508157` (0.5021039710465376` M m2^4-0.5021039710465376` M (m3^2+m23) m2^2+M^3 (0.5021039710465376` m2+0.49789602895346247` m3) m2+M^2 m1 (0.5021039710465376` m2+0.49789602895346247` m3) m2-0.49789602895346247` M (m1^2+m23) m3 m2-0.5021039710465376` M m1^2 m3^2-0.5021039710465376` M m12 (m2^2-1.` m3^2)-0.5021039710465376` m1 m12 (m2^2-1.` m3^2)-0.5021039710465376` m1 m3 (1.` m3^3-1.` m2^2 m3-1.` m23 m3-0.9916193809734177` m2 m23+m1^2 (0.9916193809734177` m2+1.` m3))) mZp^4) mZp^2+(M-m1)^2 (1/(mZp^2-8315.178393760001`)^2 0.30511302550508157` (-0.5021039710465376` m2^4+0.5021039710465376` (2 m3^2+m23) m2^2+0.995792057906925` m23 m3 m2-0.5021039710465376` m3^2 (m3^2-1.` m23)) mZp^4+(0.06066071052955202` (m2^2+2 m3 m2+m3^2-m23) (m2-m3)^2 mZp^2)/(mZp^2-8315.178393760001`)- 1.4330273763308354` (m2-m3)^2 (m2^2+2 m3 m2+m3^2-m23)) (M^2+2 m1 M+m1^2-m23)-1/(mZp^2-8315.178393760001`)^2 0.6127938468964998` ((1.` m1^2+0.5` m2^2+0.5` m3^2-1.` m12-0.4999999999999999` m23+0.9916193809734177` m2 m3) M^2+m1 (-1.` m2^2-3.9664775238936714` m3 m2-1.` m3^2+1.` m23) M+1.` m12^2-0.9999999999999998` m12 m2^2+0.5` m23^2+1.` m2^2 m3^2-0.9999999999999998` m12 m3^2-0.4999999999999999` m23 m3^2+m1^2 (0.5` m2^2+0.9916193809734178` m3 m2+0.5` m3^2-1.` m12-0.5` m23)-0.4999999999999999` m2^2 m23+1.` m12 m23-0.9916193809734177` m2 m23 m3) mZp^8) ;


Print["Proceeding with uncompiled distribution."];ChiDistributionElectrons3bodyCompiled[M_,m1_,mZp_,m12_,m23_]=ChiDistributionLeptons3body[M,m1,ElectronMass,ElectronMass,mZp,m12,m23];
ChiDistributionMuons3bodyCompiled[M_,m1_,mZp_,m12_,m23_]=ChiDistributionLeptons3body[M,m1,MuonMass,MuonMass,mZp,m12,m23];

duration=(Now-duration)//First
]


(* ::Section::Closed:: *)
(*Unitary vector combination*)


CombiningVariablesIntoAunityVectorCompiled=Null;
CompileUnityVector[]:=Module[
{duration,(*variables*)theta,z},

duration=Now;

CombiningVariablesIntoAunityVectorCompiled=
Compile[{
{theta, _Real},{z,_Real}},
{Sqrt[1-z^2]Cos[theta],Sqrt[1-z^2]Sin[theta],z},
CompilationTarget->"C",
Parallelization->True,
RuntimeAttributes->{Listable},RuntimeOptions->"Speed",
CompilationOptions->{"ExpressionOptimization" -> True,"InlineExternalDefinitions" -> True}];



(*checking if functions were compiled - if not, just return a usual mathematica function*)
If[(StringContainsQ[CompilePrint[CombiningVariablesIntoAunityVectorCompiled],"MainEvaluate"]),Print["ERROR: unity vector not compiled. Proceeding with uncompiled unity vector."];CombiningVariablesIntoAunityVectorCompiled=Null;,Print["\[Checkmark] Unitary vector compiled."];];

duration=(Now-duration)//First
]


CompileUnityVector[0]:=Null;


(* ::Chapter::Closed:: *)
(*Compiling*)


initialTimeCompilation=Now;


If[!compilationFlag,CompileBoosts[0];CompileDalitzBoundaries[0];CompileElectronPositronCollisions[0];CompileChiDecays[0];CompileUnityVector[0];,CompileBoosts[];CompileDalitzBoundaries[];CompileElectronPositronCollisions[];CompileChiDecays[];CompileUnityVector[]];


duration2=(Now-initialTimeCompilation)//First
message=StringJoin["Compilation took ", TextString[duration2], " seconds."];
temporaryMessages=CreateDialog[Style[message]];RunScheduledTask[NotebookClose[temporaryMessages];,{4,1}];


(* ::Chapter::Closed:: *)
(*Modules for random variables generation*)


(*generate random variables for 1D distributions*)
RandomVariateNew[dist_,xmin_,xmax_,nEvents_,nTrials_]:=Module[{x,ptX,ptWeights},

ptX=RandomReal[{xmin,xmax},nTrials];
ptWeights=dist[ptX];

RandomChoice[ptWeights->ptX,nEvents]
]


(*generate random variables for 2D distributions defined only in the Dalitz plane*)
RandomVariate2DDalitz[dist_,M_,m1_,m2_,m3_,nEvents_,nTrials_,Verbose_]:=Module[{m12,pt,x,y,ptXY,ptWeights,sample},

pt=ParallelTable[{x=RandomReal[{(M-m3)^2,(m1+m2)^2}],y=RandomReal[{m23minCompiled[M,m1,m2,m3,x],m23maxCompiled[M,m1,m2,m3,x]}],(m23maxCompiled[M,m1,m2,m3,x]-m23minCompiled[M,m1,m2,m3,x])*dist[x,y]},nTrials,DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}];
ptXY=Table[{pt[[i]][[1]],pt[[i]][[2]]},{i,1,nTrials}];
ptWeights=Table[pt[[i]][[3]],{i,1,nTrials}];

If[Verbose===1,sample=RandomChoice[ptWeights->ptXY,nEvents];
Print[{Show[ListPlot[ptXY,PlotStyle->Opacity[0.5,Red]]],Show[ListPlot[sample,PlotStyle->Opacity[0.8,Blue]]]}];Print[{DensityHistogram[sample],DensityPlot[dist[x,y],{x,(m1+m2)^2,(M-m3)^2},{y,(m3+m2)^2,(M-m1)^2}]}];sample,RandomChoice[ptWeights->ptXY,nEvents]]


](*specifically made for dalitz plots where we will have distributions which are functions comign with boundaries and second variable boundary depends on first variable*);


If[(StringContainsQ[CompilePrint[CombiningVariablesIntoAunityVectorCompiled],"MainEvaluate"]),ResourceFunction["RandomUnitVector"][3,10];RandomUnitVectorMine[nVectors_]:=ResourceFunction["RandomUnitVector"][3,nVectors];,RandomUnitVectorMine[nVectors_]:=CombiningVariablesIntoAunityVectorCompiled[RandomReal[{0,2\[Pi]},nVectors],RandomReal[{-1,1},nVectors]];];


(* ::Chapter::Closed:: *)
(*Modules for event generation*)


(* ::Subsection::Closed:: *)
(*Z prime generation*)


Generator["e- e+ -> \[Gamma] Z'",minPhotonAngleCM_(*in degrees*),maxPhotonAngleCM_(*in degrees*)(*photons will be generated only in this interval of angles - since the distribution grows fast at the 0 and \[Pi] the cuts speed up a lot the simulation*),velocityCM_(*3-velocity of the CM in the lab frame*),EnergyInCM_,m1_(*final state particle 1 - photon for us*),m2_(*final state particle 2 - dark photon for us*),nEvents_,nTrialsPhoton_,Verbose_]:=Module[{AmplitudeSquared,\[Theta]CMs,FourMomentaCM,FourMomentaLAB,E1,E2,p1,\[Phi]CMs,px,py,pz,vx,vy,vz,i,duration,FinalStateParticles,returning},

(*defining particles distribution from e- e+ collisions*)
AmplitudeSquared[x_]:=ElectronPositronCollisionDistributionNoEpsCompiled[EnergyInCM^2,x,m2](*the amplitude should be a function of the angle between the electron and the photon in the CM frame \[Theta]CM *);

(*initial definitions*)
FinalStateParticles={22,101};
vx=velocityCM[[1]];
vy=velocityCM[[2]];
vz=velocityCM[[3]];


(*generating the \[Theta]CM following the amplitude squared and uniform distribution for \[Phi] around the colision line*)
\[Theta]CMs=RandomVariateNew[AmplitudeSquared,minPhotonAngleCM*(\[Pi]/180),maxPhotonAngleCM*(\[Pi]/180),nEvents,nTrialsPhoton];
\[Phi]CMs=RandomReal[{0,2\[Pi]},nEvents];

If[Verbose==1,Print["Printing info about the photon angular distributions in the CM of the e+ e- collider"];Print[{Histogram[\[Theta]CMs],Plot[AmplitudeSquared[x],{x,minPhotonAngleCM*(\[Pi]/180),maxPhotonAngleCM*(\[Pi]/180)}],Histogram[\[Phi]CMs]}];,0;];


(*2 to 2 kinematics*)
E1=(EnergyInCM^2+m1^2-m2^2)/(2EnergyInCM);E2=(EnergyInCM^2+m2^2-m1^2)/(2EnergyInCM);p1=Sqrt[E1^2-m1^2];


(*assigning four momenta   -  structure {E1,p1x,p1y,p1z,E2,p2x,p2y,p2z,vx,vy,vz}*)
FourMomentaCM=AssigningMomentaElectronPositronCollisions[E1,p1,E2,vx,vy,vz,\[Phi]CMs,\[Theta]CMs,nEvents]//Transpose;

(*Boosting to the lab frame which has a velocity = -velocityCM*)
FourMomentaLAB=DoubleFourMomentaBoost@@FourMomentaCM;

ReturningParticlesElectronPositronCollisions[FinalStateParticles[[1]],FinalStateParticles[[2]],FourMomentaLAB,nEvents]
];


(* ::Subsection::Closed:: *)
(*Z prime decay (TB compiled to make it faster)*)


Generator["Z'->",ZprimeEvents_,nEvents_,mZp_,FinalStateMassesPerChannel_,ZpTotalDecayWidth_,ZpBRs_(*lighest states to heaviest : {--, -+, ++}*)]:=Module[{decayVertex,directionsCM,FourMomentaCM,FourMomentaLAB,Channel,ChannelPerDecay,E1,E2,p,p1,px,py,pz,i,k,FinalStateParticles,properLifetimes,ZpVelocity,ZprimeEvent},


(*initial definitions*)
FinalStateParticles={{"\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)"},{"\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"},{"\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)","\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)"}};


(*generating the lenght of propagation for each Zprime and decay vertex - if at the dilated lifetime of Zp times the speed of light is smaller than 0.1 cm, we consider the event as a prompt decay*)
properLifetimes=RandomVariate[ExponentialDistribution[ZpTotalDecayWidth],nEvents];
decayVertex=Table[0.197*10^-15(*position of the decay vertex in meters*)*properLifetimes[[i]](*proper lifetime*)*(ZprimeEvents[[i]][[2;;4]]/mZp)(*Zp 3-velocity times gamma which is responsible for time dilation*),{i,1,nEvents}];


(*choosing the decay channel*) (********************************************** to be Compiled*)
Channel[x_]:=Piecewise[{{1,x<=ZpBRs[[1]]},{2,ZpBRs[[1]]<x<=ZpBRs[[1]]+ZpBRs[[2]]},{3,ZpBRs[[1]]+ZpBRs[[2]]<x}}];
ChannelPerDecay=ParallelTable[Channel[RandomReal[]],{i,1,nEvents},DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}];


(*1 to 2 kinematics*)
E1[k_]:=1/(2mZp) (mZp^2+FinalStateMassesPerChannel[[k]][[1]]^2-FinalStateMassesPerChannel[[k]][[2]]^2);
E2[k_]:=1/(2mZp) (mZp^2+FinalStateMassesPerChannel[[k]][[2]]^2-FinalStateMassesPerChannel[[k]][[1]]^2);p1[k_]:=Sqrt[E1[k]^2-FinalStateMassesPerChannel[[k]][[1]]^2];


(*generating directions of the final state particle 1 - uniformly distributed in the rest-frame of Zp*)
directionsCM=RandomUnitVectorMine[nEvents];


(*assigning four momenta   -  structure {E1,p1x,p1y,p1z,E2,p2x,p2y,p2z,vx,vy,vz}*)
FourMomentaCM=(Table[{E1[ChannelPerDecay[[i]]],px=(p=(p1[ChannelPerDecay[[i]]]*directionsCM[[i]]))[[1]],py=p[[2]],pz=p[[3]],E2[ChannelPerDecay[[i]]],-px,-py,-pz,-(ZpVelocity=((ZprimeEvent=ZprimeEvents[[i]])[[2;;4]]/ZprimeEvent[[1]]))[[1]],-ZpVelocity[[2]],-ZpVelocity[[3]]},{i,1,nEvents}]//Transpose);

(*Boosting to the lab frame which has a velocity = -velocityCM*)
FourMomentaLAB=DoubleFourMomentaBoost@@FourMomentaCM;

Table[{{{FinalStateParticles[[ChannelPerDecay[[i]]]][[1]],FourMomentaLAB[[i]][[1]]},{FinalStateParticles[[ChannelPerDecay[[i]]]][[2]],FourMomentaLAB[[i]][[2]]}},decayVertex[[i]]},{i,1,nEvents}]
]


(* ::Subsection::Closed:: *)
(*SubPlus[\[Chi]] decay (TB compiled to make it faster)*)


(* ::Subsubsection::Closed:: *)
(*Generator for many decays in CM*)


(*Decays will be only into SM*)
Generator["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\) -> \!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\) \!\(\*SuperscriptBox[\(l\), \(+\)]\) \!\(\*SuperscriptBox[\(l\), \(-\)]\)",nEvents_,nTrialsChis_,mZp_,mChiPlus_,mChiMinus_,ChiPlusTotalDecayWidth_,ChiPlusBRs_(*lighest states to heaviest : {--, -+, ++}*),Verbose_]:=Module[{decayLenght,directions1CM,directions3CM,angleInCone,m1,m2,m3,m12m23PTS,E,p,E1,E3,p1abs,p2abs,p3abs,p1,p3,cos\[Theta]13,sin\[Theta]13,i,k,FinalStateParticles,vChiPlus,duration,x,y,z,\[Theta],E1i,E3i,testFunct,RotBack,d,ChiPlusDecayDistribution,Events,Channel,ChannelPerDecay,EventsNonElectrons,nEventsElectrons,nEventsMuons,nEventsPions},


(*initial definitions*)
m1=mChiMinus;


(*choosing how many events in each channel we will generate*)
Channel[x_]:=Piecewise[{{2,x<=ChiPlusBRs[[2]]}(*,{1,ChiPlusBRs[[2]]<x}}]*),{3,ChiPlusBRs[[2]]<x<=ChiPlusBRs[[2]]+ChiPlusBRs[[3]]},{1,ChiPlusBRs[[2]]+ChiPlusBRs[[3]]<x}}](*because the electron BR is the largest so may be larger than 1*);
ChannelPerDecay=ParallelTable[Channel[RandomReal[]],{i,1,nEvents},DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}];
EventsNonElectrons=Select[ChannelPerDecay,#>1.5&];
nEventsElectrons=nEvents-(EventsNonElectrons//Length);
If[nEventsElectrons<nEvents,nEventsMuons=Select[EventsNonElectrons,#<2.5&]//Length;,nEventsMuons=0;];
If[nEventsElectrons+nEventsMuons<nEvents,nEventsPions=Select[EventsNonElectrons,#>2.5&]//Length;,nEventsPions=0;];



(*Generating (decay lenghts) / ((particle velocity/c)gamma)  in the CM, so proper decay lengths*)
decayLenght=0.197*10^-15(*meters*)RandomVariate[ExponentialDistribution[ChiPlusTotalDecayWidth],nEvents];

(*    generating directions of final state particles   *)
directions1CM=RandomUnitVectorMine[nEvents](*generating directions of the final state particle 1 - uniformly distributed in the rest-frame of Chi+*);
angleInCone=RandomReal[{0,2\[Pi]},nEvents](*generating azimutal angle of particle 3 in the cone of fixed angle Subscript[\[Theta], 13] around the direction of particle 1 - uniformly distributed in the rest-frame of Chi+*);
(*for rotation things back to the right positition same as RotationMatrix[{{0,0,1},{x,y,z}}] *)
RotBack[x_,y_,z_]={{(y^2+(x^2 z)/Sqrt[x^2+y^2+z^2])/(x^2+y^2),(x y (-1+z/Sqrt[x^2+y^2+z^2]))/(x^2+y^2),x/Sqrt[x^2+y^2+z^2]},{(x y (-1+z/Sqrt[x^2+y^2+z^2]))/(x^2+y^2),(x^2+(y^2 z)/Sqrt[x^2+y^2+z^2])/(x^2+y^2),y/Sqrt[x^2+y^2+z^2]},{-(x/Sqrt[x^2+y^2+z^2]),-(y/Sqrt[x^2+y^2+z^2]),z/Sqrt[x^2+y^2+z^2]}};


(* ******************************* generating electron events ****************************************************************)

(*initial definitions*)
FinalStateParticles={"\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SuperscriptBox[\(e\), \(+\)]\)","\!\(\*SuperscriptBox[\(e\), \(-\)]\)"}; 
m2=ElectronMass;
m3=ElectronMass;


(*defining particles distribution for chi+ decays*)
ChiPlusDecayDistribution[x_,y_]:=ChiDistributionElectrons3bodyCompiled[mChiPlus,mChiMinus,mZp,x,y];


(*generating E1 and E2 in CM following the amplitude squared - using Dalitz variables m12 and m23 which are the invariant mass of particles i and j squared*) 
If[Verbose===1,duration=AbsoluteTiming[Print["Printing info about the m12 and m23 squared random values generated in the CM of the mother particle"];m12m23PTS=RandomVariate2DDalitz[ChiPlusDecayDistribution,mChiPlus,m1,m2,m3,nEventsElectrons,nTrialsChis,1];testFunct[x_,y_]:=1;Print["Printing uniform distribution to check if enough points to cover the Dalitz plane"];RandomVariate2DDalitz[testFunct,mChiPlus,m1,m2,m3,IntegerPart[2/3 nTrialsChis],nTrialsChis,1]][[1]];Print[duration];,m12m23PTS=RandomVariate2DDalitz[ChiPlusDecayDistribution,mChiPlus,m1,m2,m3,nEventsElectrons,nTrialsChis,0];](*generating decays*);


(*converting from Dalitz variables to energies*)
E1=Table[(mChiPlus^2+m1^2-m12m23PTS[[i]][[2]])/(2mChiPlus),{i,1,nEventsElectrons}];
E3=Table[(mChiPlus^2+m3^2-m12m23PTS[[i]][[1]])/(2mChiPlus),{i,1,nEventsElectrons}];



(*kinematics for 1 to 3 in CM*) 
p1abs[E_]=Sqrt[E^2-m1^2];
p2abs[E_]=Sqrt[E^2-m2^2];
p3abs[E_]=Sqrt[E^2-m3^2];
cos\[Theta]13[x_(*E1*),y_(*E3*)]=((p2abs[mChiPlus-x-y])^2-p1abs[x]^2-p3abs[y]^2)/(2p1abs[x] p3abs[y]);
sin\[Theta]13[x_(*E1*),y_(*E3*)]=\[Sqrt](1-(((p2abs[mChiPlus-x-y])^2-p1abs[x]^2-p3abs[y]^2)/(2p1abs[x] p3abs[y]))^2)(*since Subscript[\[Theta], 13] from 0 to \[Pi], the sine is always positive*);



directions3CM=Table[RotBack[(d=directions1CM[[i]])[[1]],d[[2]],d[[3]]] . {Cos[(\[Theta]=angleInCone[[i]])]*sin\[Theta]13[(E1i=E1[[i]]),(E3i=E3[[i]])],Sin[\[Theta]]*sin\[Theta]13[E1i,E3i],cos\[Theta]13[E1i,E3i]},{i,1,nEventsElectrons}](*generating directions of the final state particle 3*);


(*    creating the four-momenta of the final states   *)
Events=Table[{{{FinalStateParticles[[1]],Flatten[{E1[[i]],p1=(p1abs[E1[[i]]]*directions1CM[[i]])}]},{FinalStateParticles[[3]],Flatten[{E3[[i]],p3=(p3abs[E3[[i]]]*directions3CM[[i]])}]},{FinalStateParticles[[2]],Flatten[{mChiPlus-E1[[i]]-E3[[i]],-(p1+p3)}]}},decayLenght[[i]]},{i,1,nEventsElectrons}](*assigning four momenta in the rest-frame of chi+*);


(* ******************************* generating muon events ****************************************************************)
If[nEventsMuons===0,Return[Events,Module],0;];

(*initial definitions*)
FinalStateParticles={"\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SuperscriptBox[\(\[Mu]\), \(+\)]\)","\!\(\*SuperscriptBox[\(\[Mu]\), \(-\)]\)"}; 
m2=MuonMass;
m3=MuonMass;


(*defining particles distribution for chi+ decays*)
ChiPlusDecayDistribution[x_,y_]:=ChiDistributionMuons3bodyCompiled[mChiPlus,mChiMinus,mZp,x,y];


(*generating E1 and E2 in CM following the amplitude squared - using Dalitz variables m12 and m23 which are the invariant mass of particles i and j squared*) 
If[Verbose===1,duration=AbsoluteTiming[Print["Printing info about the m12 and m23 squared random values generated in the CM of the mother particle"];m12m23PTS=RandomVariate2DDalitz[ChiPlusDecayDistribution,mChiPlus,m1,m2,m3,nEventsMuons,nTrialsChis,1];testFunct[x_,y_]:=1;Print["Printing uniform distribution to check if enough points to cover the Dalitz plane"];RandomVariate2DDalitz[testFunct,mChiPlus,m1,m2,m3,IntegerPart[2/3 nTrialsChis],nTrialsChis,1]][[1]];Print[duration];,m12m23PTS=RandomVariate2DDalitz[ChiPlusDecayDistribution,mChiPlus,m1,m2,m3,nEventsMuons,nTrialsChis,0];](*generating decays*);


(*converting from Dalitz variables to energies*)
E1=Table[(mChiPlus^2+m1^2-m12m23PTS[[i]][[2]])/(2mChiPlus),{i,1,nEventsMuons}];
E3=Table[(mChiPlus^2+m3^2-m12m23PTS[[i]][[1]])/(2mChiPlus),{i,1,nEventsMuons}];



(*kinematics for 1 to 3 in CM*) 
p1abs[E_]=Sqrt[E^2-m1^2];
p2abs[E_]=Sqrt[E^2-m2^2];
p3abs[E_]=Sqrt[E^2-m3^2];
cos\[Theta]13[x_(*E1*),y_(*E3*)]=((p2abs[mChiPlus-x-y])^2-p1abs[x]^2-p3abs[y]^2)/(2p1abs[x] p3abs[y]);
sin\[Theta]13[x_(*E1*),y_(*E3*)]=\[Sqrt](1-(((p2abs[mChiPlus-x-y])^2-p1abs[x]^2-p3abs[y]^2)/(2p1abs[x] p3abs[y]))^2)(*since Subscript[\[Theta], 13] from 0 to \[Pi], the sine is always positive*);



directions3CM=Table[RotBack[(d=directions1CM[[i+nEventsElectrons]])[[1]],d[[2]],d[[3]]] . {Cos[(\[Theta]=angleInCone[[i+nEventsElectrons]])]*sin\[Theta]13[(E1i=E1[[i]]),(E3i=E3[[i]])],Sin[\[Theta]]*sin\[Theta]13[E1i,E3i],cos\[Theta]13[E1i,E3i]},{i,1,nEventsMuons}](*generating directions of the final state particle 3*);


(*    creating the four-momenta of the final states   *)
Events=Join[Events,Table[{{{FinalStateParticles[[1]],Flatten[{E1[[i]],p1=(p1abs[E1[[i]]]*directions1CM[[i+nEventsElectrons]])}]},{FinalStateParticles[[3]],Flatten[{E3[[i]],p3=(p3abs[E3[[i]]]*directions3CM[[i]])}]},{FinalStateParticles[[2]],Flatten[{mChiPlus-E1[[i]]-E3[[i]],-(p1+p3)}]}},decayLenght[[i+nEventsElectrons]]},{i,1,nEventsMuons}]];(*assigning four momenta in the rest-frame of chi+*)

(* ******************************* generating pions events ****************************************************************)
If[nEventsPions===0,Return[Events,Module],0;];

(*initial definitions*)
FinalStateParticles={"\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)","\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)","\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)"}; 
m2=PionMass;
m3=PionMass;


(*defining particles distribution for chi+ decays*)
ChiPlusDecayDistribution[x_,y_]:=ChiDistributionPions3bodyCompiled[mChiPlus,mChiMinus,mZp,x,y];


(*generating E1 and E2 in CM following the amplitude squared - using Dalitz variables m12 and m23 which are the invariant mass of particles i and j squared*) 
If[Verbose===1,duration=AbsoluteTiming[Print["Printing info about the m12 and m23 squared random values generated in the CM of the mother particle"];m12m23PTS=RandomVariate2DDalitz[ChiPlusDecayDistribution,mChiPlus,m1,m2,m3,nEventsPions,nTrialsChis,1];testFunct[x_,y_]:=1;Print["Printing uniform distribution to check if enough points to cover the Dalitz plane"];RandomVariate2DDalitz[testFunct,mChiPlus,m1,m2,m3,IntegerPart[2/3 nTrialsChis],nTrialsChis,1]][[1]];Print[duration];,m12m23PTS=RandomVariate2DDalitz[ChiPlusDecayDistribution,mChiPlus,m1,m2,m3,nEventsPions,nTrialsChis,0];](*generating decays*);


(*converting from Dalitz variables to energies*)
E1=Table[(mChiPlus^2+m1^2-m12m23PTS[[i]][[2]])/(2mChiPlus),{i,1,nEventsPions}];
E3=Table[(mChiPlus^2+m3^2-m12m23PTS[[i]][[1]])/(2mChiPlus),{i,1,nEventsPions}];



(*kinematics for 1 to 3 in CM*) 
p1abs[E_]=Sqrt[E^2-m1^2];
p2abs[E_]=Sqrt[E^2-m2^2];
p3abs[E_]=Sqrt[E^2-m3^2];
cos\[Theta]13[x_(*E1*),y_(*E3*)]=((p2abs[mChiPlus-x-y])^2-p1abs[x]^2-p3abs[y]^2)/(2p1abs[x] p3abs[y]);
sin\[Theta]13[x_(*E1*),y_(*E3*)]=Sqrt[1-(((p2abs[mChiPlus-x-y])^2-p1abs[x]^2-p3abs[y]^2)/(2p1abs[x] p3abs[y]))^2](*since Subscript[\[Theta], 13] from 0 to \[Pi], the sine is always positive*);



directions3CM=Table[RotBack[(d=directions1CM[[i+nEventsElectrons+nEventsMuons]])[[1]],d[[2]],d[[3]]] . {Cos[(\[Theta]=angleInCone[[i+nEventsElectrons+nEventsMuons]])]*sin\[Theta]13[(E1i=E1[[i]]),(E3i=E3[[i]])],Sin[\[Theta]]*sin\[Theta]13[E1i,E3i],cos\[Theta]13[E1i,E3i]},{i,1,nEventsPions}](*generating directions of the final state particle 3*);


(*    creating the four-momenta of the final states   *)
Join[Events,Table[{{{FinalStateParticles[[1]],Flatten[{E1[[i]],p1=(p1abs[E1[[i]]]*directions1CM[[i+nEventsElectrons+nEventsMuons]])}]},{FinalStateParticles[[3]],Flatten[{E3[[i]],p3=(p3abs[E3[[i]]]*directions3CM[[i]])}]},{FinalStateParticles[[2]],Flatten[{mChiPlus-E1[[i]]-E3[[i]],-(p1+p3)}]}},decayLenght[[i+nEventsElectrons+nEventsMuons]]},{i,1,nEventsPions}]](*assigning four momenta in the rest-frame of chi+*)
]


(* ::Subsubsection::Closed:: *)
(*Boosting to LAB for single decay*)


SingleEventBoost["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\) -> \!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\) \!\(\*SuperscriptBox[\(l\), \(+\)]\) \!\(\*SuperscriptBox[\(l\), \(-\)]\)",ChiPlusEvent_,DecayEventCM_,mChiPlus_]:=Module[{decayVertex,pChiPlus,vChiPlus,BoostedFourMomenta},
pChiPlus=ChiPlusEvent[[1]][[2]][[2;;4]](*chi+ momentum*);
vChiPlus=(pChiPlus/ChiPlusEvent[[1]][[2]][[1]])(*chi+ velocity*);

(*generating the lenght of propagation for each ChiPlus and decay vertex - if the dilated lifetime of Chi+ times the speed of light is smaller than 0.1 cm, we consider the event as a prompt decay*)
decayVertex=DecayEventCM[[2]](*proper decay lenght/velocity*)*(pChiPlus/mChiPlus)(*Zp 3-velocity times gamma*)+ChiPlusEvent[[2]];


(*Boosting to the lab frame which has a velocity = - ChiPlusVelocity*)
BoostedFourMomenta=TripleFourMomentaBoost[DecayEventCM[[1]][[1]][[2]][[1]],DecayEventCM[[1]][[1]][[2]][[2]],DecayEventCM[[1]][[1]][[2]][[3]],DecayEventCM[[1]][[1]][[2]][[4]],DecayEventCM[[1]][[2]][[2]][[1]],DecayEventCM[[1]][[2]][[2]][[2]],DecayEventCM[[1]][[2]][[2]][[3]],DecayEventCM[[1]][[2]][[2]][[4]],DecayEventCM[[1]][[3]][[2]][[1]],DecayEventCM[[1]][[3]][[2]][[2]],DecayEventCM[[1]][[3]][[2]][[3]],DecayEventCM[[1]][[3]][[2]][[4]],-vChiPlus[[1]],-vChiPlus[[2]],-vChiPlus[[3]]];

(*returning a list with the form {{{particle name,its four momenta},{particle name,its four momenta},{particle name,its four momenta}},Chi+ decay vertex position}*)
{{{DecayEventCM[[1]][[1]][[1]],BoostedFourMomenta[[1]]},{DecayEventCM[[1]][[2]][[1]],BoostedFourMomenta[[2]]},{DecayEventCM[[1]][[3]][[1]],BoostedFourMomenta[[3]]}},decayVertex}
];


(* ::Subsection::Closed:: *)
(*Full event generator (TB compiled to make it faster)*)


(* ::Subsubsection::Closed:: *)
(*Pre definitions*)


DecayZpDaugthers[ZpDecayEvent_,{DecayEvent1CM_,DecayEvent2CM_},mChiPlus_]:=Module[{Particle1,Particle2,m2,m3,m23min,m23max,ElectronPositronCollisionDistributionForMax,ElectronPositronCollisionDistributionMax},


Particle1=If[ZpDecayEvent[[1]][[1]][[1]]==="\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",SingleEventBoost["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\) -> \!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\) \!\(\*SuperscriptBox[\(l\), \(+\)]\) \!\(\*SuperscriptBox[\(l\), \(-\)]\)",{ZpDecayEvent[[1]][[1]],ZpDecayEvent[[2]]},DecayEvent1CM,mChiPlus],Null];
Particle2=If[ZpDecayEvent[[1]][[2]][[1]]==="\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",SingleEventBoost["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\) -> \!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\) \!\(\*SuperscriptBox[\(l\), \(+\)]\) \!\(\*SuperscriptBox[\(l\), \(-\)]\)",{ZpDecayEvent[[1]][[2]],ZpDecayEvent[[2]]},DecayEvent2CM,mChiPlus],Null];

{Particle1,Particle2}
]

EventCombiner[CollisionEvent_,ZpDecay_,ChiPlusDecay_]:=Module[{Photon,Zp,ZpDecayVertex,ZpDaugther1,ZpDaugther2,FullEvent},

Photon=CollisionEvent[[1]];
Zp=CollisionEvent[[2]];
ZpDecayVertex=ZpDecay[[2]];
ZpDaugther1=ZpDecay[[1]][[1]];
ZpDaugther2=ZpDecay[[1]][[2]];


If[ChiPlusDecay[[1]]===Null&&ChiPlusDecay[[2]]===Null,FullEvent={Photon,{Zp,{ZpDaugther1,ZpDaugther2}}},0;];If[ChiPlusDecay[[1]]===Null&&(ChiPlusDecay[[2]]//Length)>1,FullEvent={Photon,{Zp,{{ZpDaugther1,{ZpDaugther2,ChiPlusDecay[[2]]}},ZpDecayVertex}}},0;];If[(ChiPlusDecay[[1]]//Length)>1&&(ChiPlusDecay[[2]]//Length)>1,FullEvent={Photon,{Zp,{{{ZpDaugther1,ChiPlusDecay[[1]]},{ZpDaugther2,ChiPlusDecay[[2]]}},ZpDecayVertex}}},0;];
FullEvent
]


(* ::Subsubsection::Closed:: *)
(*Main definition*)


(*parameters for generation should come in a list of the type {mDM,mZp,\[Alpha]D,\[Epsilon],\[Delta]y,\[CapitalDelta]f} - all masses in GeV*)
FullEventGenerator["e- e+ -> \[Gamma] Z'",ParametersList_,minPhotonAngleCM_(*in degrees*),maxPhotonAngleCM_(*in degrees*),EnergyCoM_(*in GeV*),velocityCoM_(*normalized by the speed of light*),nEvents_,nTrialsPhoton_,nTrialsChi_,verbose_]:=Module[{x,y,m12,m23,mDM,m1(*mDM so not actually used*),m23min,m23max,mZp,mChiPlus,ZpTotalDecayWidth,ZpBranchingRatios,ChiPlusTotalDecayWidth,ChiPlusBranchingRatios,BRelectrons,ElectronPositronCollisionDistribution,ChiPlusDecayDistribution,ChiPlusDecayDistributionForMax,ChiPlusDecayDistributionMax,distributionPlots,(*event related variables :*)CollisionEvents,ZpDecays,ChiPlusDecaysCM,ChiPlusDecays,FullEvents,(*to be removed*)CollisionEventstemp},

(*defining variables related to the model*)
mDM=ParametersList[[1]];mZp=ParametersList[[2]];mChiPlus=(1+ParametersList[[6]])*mDM;


(*computing particles decay widths*)
ZpTotalDecayWidth=\[CapitalGamma]["Z'",ParametersList];
ZpBranchingRatios={\[CapitalGamma]["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)",ParametersList]/ZpTotalDecayWidth,\[CapitalGamma]["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",ParametersList]/ZpTotalDecayWidth,\[CapitalGamma]["Z'->\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",ParametersList]/ZpTotalDecayWidth};
ChiPlusTotalDecayWidth=Mean[ParallelTable[\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",ParametersList],4,DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}]];
ChiPlusBranchingRatios=If[ParametersList[[6]]*mDM<2*MuonMass,{1,0,0},{Mean[ParallelTable[\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(e\), \(-\)]\)\!\(\*SuperscriptBox[\(e\), \(+\)]\)",ParametersList],11,DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}]]/ChiPlusTotalDecayWidth,Mean[ParallelTable[\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Mu]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Mu]\), \(+\)]\)",ParametersList],11,DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}]]/ChiPlusTotalDecayWidth,Mean[ParallelTable[\[CapitalGamma]["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)->\!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(-\)]\)\!\(\*SuperscriptBox[\(\[Pi]\), \(+\)]\)",ParametersList],11,DistributedContexts->{"eeCOLLIDERSeventGENERATORniDMv3`Private`"}]]/ChiPlusTotalDecayWidth}];(*14 and 11 for the average is chosen in order to have less than 1% error in the branching ratio for 10 times trials of 1000 computations each trial*)  
If[!pionsFlag,ChiPlusBranchingRatios={ChiPlusBranchingRatios[[1]]/(ChiPlusBranchingRatios[[1]]+ChiPlusBranchingRatios[[2]]),ChiPlusBranchingRatios[[2]]/(ChiPlusBranchingRatios[[1]]+ChiPlusBranchingRatios[[2]]),0};,0;];

If[verbose===1,Print["\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(Z'\)]\) = ",ZpTotalDecayWidth, " with \!\(\*SubscriptBox[\(Br\), \(Z'\)]\) = ", ZpBranchingRatios, "; \!\(\*SubscriptBox[\(\[CapitalGamma]\), \(\[Chi]\)]\) = ",ChiPlusTotalDecayWidth, " with \!\(\*SubscriptBox[\(Br\), \(\[Chi]\)]\) = ", ChiPlusBranchingRatios];,0;];



(*start generating events ********************************************************************************************)
CollisionEventstemp=Generator["e- e+ -> \[Gamma] Z'",minPhotonAngleCM,maxPhotonAngleCM,velocityCoM,EnergyCoM,0(*photon mass*),mZp,nEvents,nTrialsPhoton,verbose];
(*to be removed once everything is compiled and so we dont use names as strings anymore*)CollisionEvents=Table[{{"\[Gamma]",CollisionEventstemp[[i]][[1]][[2]]},{"Z'",CollisionEventstemp[[i]][[2]][[2]]}},{i,1,nEvents}];


ZpDecays=Generator["Z'->",Table[CollisionEvents[[i]][[2]][[2]],{i,1,nEvents}],nEvents,mZp,{{mDM,mDM},{mDM,mChiPlus},{mChiPlus,mChiPlus}},ZpTotalDecayWidth,ZpBranchingRatios];

(*info about chi+ produced*)
If[verbose===1,Print[Histogram[Table[If[ZpDecays[[i]][[1]][[2]][[1]]==="\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\)",ZpDecays[[i]][[1]][[2]][[2]][[1]],Null],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(+\)]\)"]];,0;](*info about chi+ since no verbose for it - plotting only energies of one of ch+ decays per event*);


ChiPlusDecaysCM=Partition[RandomSample[Generator["\!\(\*SubscriptBox[\(\[Chi]\), \(+\)]\) -> \!\(\*SubscriptBox[\(\[Chi]\), \(-\)]\) \!\(\*SuperscriptBox[\(l\), \(+\)]\) \!\(\*SuperscriptBox[\(l\), \(-\)]\)",2(*number of chis produced per Zp decay*)*nEvents,nTrialsChi,mZp,mChiPlus,mDM,ChiPlusTotalDecayWidth,ChiPlusBranchingRatios,verbose]],2];

ChiPlusDecays=Table[DecayZpDaugthers[ZpDecays[[i]],ChiPlusDecaysCM[[i]],mChiPlus],{i,1,nEvents}];

If[verbose===1,Print[{DensityHistogram[Table[If[(ChiPlusDecays[[i]][[2]]//Length)>1,{ChiPlusDecays[[i]][[2]][[1]][[1]][[2]][[1]],ChiPlusDecays[[i]][[2]][[1]][[3]][[2]][[1]]},Null],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(1\)]\) x \!\(\*SubscriptBox[\(E\), \(3\)]\)"],Histogram[Table[If[(ChiPlusDecays[[i]][[2]]//Length)>1,ChiPlusDecays[[i]][[2]][[1]][[2]][[2]][[1]],Null],{i,1,nEvents}],PlotLabel->"\!\(\*SubscriptBox[\(E\), \(2\)]\)"]}];,0;](*info about chi+ decays since verbose off for it - plotting only energies of one of ch+ decays per event*);



FullEvents=Table[EventCombiner[CollisionEvents[[i]],ZpDecays[[i]],ChiPlusDecays[[i]]],{i,1,nEvents}](*all events combined to be returned*);

If[verbose===1,Print[FullEventsDisplay["e- e+ -> \[Gamma] Z'",FullEvents[[1]],1,1]];,0;](*show an event to the user*);

FullEvents
]


(* ::Subsection::Closed:: *)
(*Generating only SM particles (TB compiled to make it faster)*)


(* ::Text:: *)
(*TBC atomic problem - not fully relevant but find out origin*)


OnlySMEventSelected["e- e+ -> \[Gamma] Z'",EventsWithFullInfo_,verbose_]:=Module[{nEvents,Photons,Vertex1NegativeChargedParticle,Vertex1PositiveChargedParticle,Vertex1Position,Vertex2NegativeChargedParticle,Vertex2PositiveChargedParticle,Vertex2Position,j,jmax},



nEvents=EventsWithFullInfo//Length;
If[verbose===1,Print[nEvents];,0;];

Photons=Table[EventsWithFullInfo[[i]][[1]],{i,1,nEvents}];


(*data from vertex 1*)
Vertex1Position=Table[If[(EventsWithFullInfo[[i]][[2]][[2]][[1]][[1]]//Flatten//Length)>5,EventsWithFullInfo[[i]][[2]][[2]][[1]][[1]][[2]][[2]],Null],{i,1,nEvents}];
Vertex1NegativeChargedParticle=Table[If[(EventsWithFullInfo[[i]][[2]][[2]][[1]][[1]]//Flatten//Length)>5,EventsWithFullInfo[[i]][[2]][[2]][[1]][[1]][[2]][[1]][[2]],Null],{i,1,nEvents}];
Vertex1PositiveChargedParticle=Table[If[(EventsWithFullInfo[[i]][[2]][[2]][[1]][[1]]//Flatten//Length)>5,EventsWithFullInfo[[i]][[2]][[2]][[1]][[1]][[2]][[1]][[3]],Null],{i,1,nEvents}];

(*data from vertex 2*)
Vertex2Position=Table[If[(EventsWithFullInfo[[i]][[2]][[2]][[1]][[2]]//Flatten//Length)>5,EventsWithFullInfo[[i]][[2]][[2]][[1]][[2]][[2]][[2]],Null],{i,1,nEvents}];
Vertex2NegativeChargedParticle=Table[If[(EventsWithFullInfo[[i]][[2]][[2]][[1]][[2]]//Flatten//Length)>5,EventsWithFullInfo[[i]][[2]][[2]][[1]][[2]][[2]][[1]][[2]],Null],{i,1,nEvents}];
Vertex2PositiveChargedParticle=Table[If[(EventsWithFullInfo[[i]][[2]][[2]][[1]][[2]]//Flatten//Length)>5,EventsWithFullInfo[[i]][[2]][[2]][[1]][[2]][[2]][[1]][[3]],Null],{i,1,nEvents}];

(*info about SM particles produced*)
If[verbose===1,Print[Histogram[Table[Photons[[i]][[2]][[1]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(\[Gamma]\)]\)"]];Print[{Histogram[Table[If[Vertex1Position[[i]]===Null,Null,Norm[Vertex1Position[[i]]]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(\[Chi]\), \(+1\)]\) propagation lenght"],Histogram[Table[If[Vertex1NegativeChargedParticle[[i]]===Null,Null,Vertex1NegativeChargedParticle[[i]][[2]][[1]]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(\(1\)\(-\)\)]\)"],Histogram[Table[If[Vertex1PositiveChargedParticle[[i]]===Null,Null,Vertex1PositiveChargedParticle[[i]][[2]][[1]]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(\(1\)\(+\)\)]\)"]}];Print[{Histogram[Table[If[Vertex2Position[[i]]===Null,Null,Norm[Vertex2Position[[i]]]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(\[Chi]\), \(+1\)]\) propagation lenght"],Histogram[Table[If[Vertex2NegativeChargedParticle[[i]]===Null,Null,Vertex2NegativeChargedParticle[[i]][[2]][[1]]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(\(1\)\(-\)\)]\)"],Histogram[Table[If[Vertex2PositiveChargedParticle[[i]]===Null,Null,Vertex2PositiveChargedParticle[[i]][[2]][[1]]],{i,1,nEvents}],ImageSize->Small,PlotLabel->"\!\(\*SubscriptBox[\(E\), \(\(1\)\(+\)\)]\)"]}];,0;];



Table[{Photons[[i]],If[Vertex1Position[[i]]===Null,Null,{Vertex1Position[[i]],Vertex1NegativeChargedParticle[[i]],Vertex1PositiveChargedParticle[[i]]}],If[Vertex2Position[[i]]===Null,Null,{Vertex2Position[[i]],Vertex2NegativeChargedParticle[[i]],Vertex2PositiveChargedParticle[[i]]}]},{i,1,nEvents}]
]


(* ::Title::Closed:: *)
(*End*)


duration=(Now-initialTime)//First


Print["Time time elapsed for loading the package was ", duration," seconds.\n\n If no error in compilation happened, the following example should run in about 1-5 seconds (from third run and so on) for a computer with 4 cores available (check your number of cores with $ProcessorCount):\n\t\tFullEventGenerator[\"e- e+ -> \[Gamma] Z'\",{0.25,3*0.25,0.1,3*\!\(\*SuperscriptBox[\(10\), \(-3\)]\),0,0.4},32.5,99,10.5641,{0,0,0.487603},3000,20000,20000,0];//AbsoluteTiming \n\nRun the code once with full verbose to check if the results seem reasonable:\n\t\tSMEventsDisplay[\"e- e+ -> \[Gamma] Z'\",OnlySMEventSelected[\"e- e+ -> \[Gamma] Z'\",FullEventGenerator[\"e- e+ -> \[Gamma] Z'\",{0.25,3*0.25,0.1,3*\!\(\*SuperscriptBox[\(10\), \(-3\)]\),0,0.4},32.5,99,10.5641,{0,0,0.487603},3000,20000,20000,1],1][[1]],1]\n\t\tFullEventsDisplay[\"e- e+ -> \[Gamma] Z'\",FullEventGenerator[\"e- e+ -> \[Gamma] Z'\",{0.25,3*0.25,0.1,3*\!\(\*SuperscriptBox[\(10\), \(-3\)]\),0,0.4},32.5,99,10.5641,{0,0,0.487603},100,10000,10000,0][[1]],1,0]\n\n"];


End[]


EndPackage[];
