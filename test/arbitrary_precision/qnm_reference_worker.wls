(* Independent Cook-Zalutskiy characteristic-exponent radial recurrence,
   evaluated with iterative modified Lentz; no QNM paclet CF calls in method2.
   Equations (31),(44), arXiv:1410.7698; angular lambda from BHPT SWSH.
   Method1 is separately implemented QNM 0.4.1 Leaver routine at precision80.
   Roots follow their explicit seeds, never the inversion number alone. *)
Needs["QNM`"];Needs["Teukolsky`"];Needs["SpinWeightedSpheroidalHarmonics`"];
idx=ToExpression[$ScriptCommandLine[[2]]];out=OpenWrite[$ScriptCommandLine[[3]],PageWidth->Infinity];
emit[k_,v_]:=(WriteString[out,tag,"|",k,"|",ToString[InputForm[v]],"\n"];Flush[out]);
all=Join[Flatten[Table[{a,n},{a,{2/5,9/10}},{n,0,9}],1],{{2/5,"8_1"},{9/10,"8_1"}}];
{a0,label}=all[[idx]];n=If[label==="8_1",9,label];
tag="a="<>ToString[InputForm[a0]]<>";l=2;m=2;n="<>ToString[label];
seeds04={.4398419217-.086881962 I,.4208466772-.264733449 I,.3877792547-.453259499 I,.3488558646-.653808114 I,.3115460811-.861738034 I,.2816910047-1.071525405 I,.2635708064-1.283473343 I,.256-1.50 I,.252-1.72 I,.245-2.16 I};
seeds09={.6716142721-.064869236 I,.66765755-.19525207 I,.65982668-.32751838 I,.64786902-.46286427 I,.62983602-.60329508 I,.53691672-.74872346 I,.60222947-.77034094 I,.62-.94 I,.6253445749-1.0819705743 I,.63-1.38 I};
seed=If[label==="8_1",If[a0==2/5,.25-1.94 I,.6282699784-1.2300085359 I],If[a0==2/5,seeds04[[n+1]],seeds09[[n+1]]]];
emit["seed",seed];emit["inversion",n];emit["label_status","seed-selected; Schwarzschild continuation not yet certified"];
$RecursionLimit=65536;
cf[ww_?NumericQ,aa0_?NumericQ,inv_Integer,tol_:10^-90]:=Module[{w=SetPrecision[ww,220],aa=N[aa0,220],s=-2,m=2,k,rp,rm,z,xi,eta,p,al,ga,de,sg,A,d0,d1,d2,d3,d4,an,bn,cn,left=0,ii,j,aj,bj,f,c,d,del,tiny=10^-230},
 k=Sqrt[1-aa^2];rp=1+k;rm=1-k;z=I w;xi=-s-I(2 w rp-m aa)/(2 k);eta=-I(2 w rm-m aa)/(2 k);
 p=k z;al=1+2s+xi+eta-2z;ga=1+s+2eta;de=1+s+2xi;
 A=SpinWeightedSpheroidalEigenvalue[s,2,m,aa w]+2 aa m w-aa^2 w^2;
 sg=A+aa^2 w^2-8w^2+p(2al+ga-de)+(1+s-(ga+de)/2)(s+(ga+de)/2);
 {d0,d1,d2,d3,d4}={de,4p-2al+ga-de-2,2al-ga+2,al(4p-de)-sg,al(al-ga+1)};
 an[q_]:=q^2+(d0+1)q+d0;bn[q_]:=-2q^2+(d1+2)q+d3;cn[q_]:=q^2+(d2-3)q+d4-d2+2;
 Do[left=an[ii]/(bn[ii]-cn[ii]left),{ii,0,inv-1}];
 f=N[tiny,220];c=f;d=0;j=0;del=0;
 While[j<200000 && (j<100 || Abs[del-1]>tol),j++;ii=inv+j;
 aj=-an[ii-1]/cn[ii-1];bj=bn[ii]/cn[ii];
 d=bj+aj d;If[d==0,d=tiny];c=bj+aj/c;If[c==0,c=tiny];
 (* Fixed internal precision prevents cumulative pessimistic significance
    tracking from treating a nonzero Lentz update as zero. Accuracy is
    accepted only by agreement with the independently computed MST root. *)
 d=SetPrecision[1/d,220];c=SetPrecision[c,220];del=SetPrecision[c d,220];f=SetPrecision[f del,220]];
 lastCF={Abs[del-1],j};If[j>=200000,Return[Indeterminate]];
 bn[inv]-cn[inv]left+cn[inv]f
];
t0=AbsoluteTime[];
res=TimeConstrained[
 w1=$Failed;emit["QNM_paclet_Leaver_status","not used: n0 diagnostic failed with line-search warnings and only 10-digit seed agreement; see qnm_smoke.txt"];
 a=a0;fun[w_?NumericQ]:=cf[w,a,n];
 emit["method2.begin",AbsoluteTime[]];
 w2=w/.FindRoot[fun[w]==0,{w,SetPrecision[seed,100]},WorkingPrecision->100,AccuracyGoal->80,PrecisionGoal->80,MaxIterations->60];
 emit["omega.method2.independent_Leaver100",w2];emit["abs_cf",Abs[fun[w2]]];emit["cf_last_change_and_iterations",lastCF];
 emit["digits.method1_method2",If[NumberQ[w1]&&NumberQ[w2],-Log[10,Abs[w2-w1]/Abs[w2]],"INVALID"]];
 emit["method1.MST.begin",AbsoluteTime[]];
 wMST=QNMFrequency[-2,2,2,n,N[a0,120],Method->{"IncidenceAmplitude","InitialGuess"->w2}];
 emit["omega.method1.QNM_MST120",wMST];
 emit["digits.MST120_independentCF100",If[NumberQ[wMST]&&NumberQ[w2],N[-Log[10,Abs[(SetPrecision[wMST,Infinity]-SetPrecision[w2,Infinity])/SetPrecision[wMST,Infinity]]],15],"INVALID"]];
 a=N[a0,120];wr=SetPrecision[wMST,120];
 fIn[ww_?NumericQ]:=Module[{r,amp},r=TeukolskyRadial[-2,2,2,a,ww,"BoundaryConditions"->"In",WorkingPrecision->120];amp=r["Amplitudes"];amp["Incidence"]/amp["Transmission"]];
 deriv[h_]:=(fIn[wr-2h]-8 fIn[wr-h]+8 fIn[wr+h]-fIn[wr+2h])/(12h);
 emit["derivative.begin",AbsoluteTime[]];
 d15=deriv[N[10^-15 Abs[wr],120]];emit["dA_inc.h1e-15",d15];
 d12=deriv[N[10^-12 Abs[wr],120]];emit["dA_inc.h1e-12",d12];
 emit["digits.derivative_steps",-Log[10,Abs[d15-d12]/Abs[d15]]];
 av=fIn[wr];emit["A_inc.at_root",av];emit["root_error_from_amplitude",Abs[av]/Abs[d15]];
 "completed",1800,$TimedOut];
emit["status",res];emit["wall_seconds",AbsoluteTime[]-t0];Close[out];Quit[];
