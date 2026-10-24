(* Independent arbitrary-precision references for the arbitrary-precision GSN acceptance tests, from the
   Black Hole Perturbation Toolkit Teukolsky paclet (MST method). For every case: the separation constant,
   the renormalized angular momentum, the In and Up amplitude ratios (normalization independent) and the
   In/Up radial functions and their r-derivatives at r' = 3.1126842514806725, all divided by the respective
   transmission amplitude. Output: one record per line, "case|key|value", values as InputForm at the stated
   precision. Usage: wolframscript -file bhpt_acceptance_reference.wls <workingprecision> <group>
   group = "a09l2" | "l3m1" | "axis" | "a04l2" *)
Needs["Teukolsky`"];
Needs["SpinWeightedSpheroidalHarmonics`"];
prec = ToExpression[$ScriptCommandLine[[2]]];
group = $ScriptCommandLine[[3]];
$RecursionLimit = 65536;
rprime = SetPrecision[3.1126842514806725, prec];
realFreqs = {1/20, 1/2, 2, 10};
complexFreqs = {1/2 - 3/10 I, 3/10 - 15/10 I, 1/2 + 3/10 I, 2 - 2 I};
axisFreqs = Flatten[Table[{10^-30 - sigma I, -10^-30 - sigma I, sigma I}, {sigma, {1/10, 1/2, 1, 3}}]];
cases = Switch[group,
  "a09l2", Table[{9/10, 2, 2, w}, {w, Join[realFreqs, complexFreqs]}],
  "l3m1", Join[Table[{4/10, 3, 1, w}, {w, Join[realFreqs, complexFreqs]}], Table[{9/10, 3, 1, w}, {w, Join[realFreqs, complexFreqs]}]],
  "axis", Join[Table[{4/10, 2, 2, w}, {w, axisFreqs}], Table[{9/10, 2, 2, w}, {w, axisFreqs}]],
  "a04l2", Table[{4/10, 2, 2, w}, {w, Join[realFreqs, complexFreqs]}],
  "omega10", {{4/10, 2, 2, 10}, {9/10, 2, 2, 10}, {4/10, 3, 1, 10}, {9/10, 3, 1, 10}, {4/10, 2, 2, 2 - 2 I}, {9/10, 2, 2, 2 - 2 I}},
  _, {}];
emit[tag_, key_, val_] := Print[tag, "|", key, "|", ToString[InputForm[val]]];
Do[
  CheckAbort[Module[{a0, l, m, w0, a, w, lam, Rin, Rup, ai, au, tag, t0},
    {a0, l, m, w0} = c;
    tag = StringJoin["a=", ToString[InputForm[a0]], ";l=", ToString[l], ";m=", ToString[m], ";omega=", ToString[InputForm[w0]]];
    t0 = AbsoluteTime[];
    a = SetPrecision[a0, prec]; w = SetPrecision[w0, prec];
    lam = SpinWeightedSpheroidalEigenvalue[-2, l, m, a w];
    emit[tag, "lambda", lam];
    Rin = Quiet[TeukolskyRadial[-2, l, m, a, w, "BoundaryConditions" -> "In", WorkingPrecision -> prec]];
    ai = Rin["Amplitudes"];
    emit[tag, "nu", Rin["RenormalizedAngularMomentum"]];
    emit[tag, "In.Incidence/Transmission", ai["Incidence"]/ai["Transmission"]];
    emit[tag, "In.Reflection/Transmission", ai["Reflection"]/ai["Transmission"]];
    emit[tag, "In.R(rprime)/Transmission", Rin[rprime]/ai["Transmission"]];
    emit[tag, "In.dR(rprime)/Transmission", Rin'[rprime]/ai["Transmission"]];
    Rup = Quiet[TeukolskyRadial[-2, l, m, a, w, "BoundaryConditions" -> "Up", WorkingPrecision -> prec]];
    au = Rup["Amplitudes"];
    emit[tag, "Up.Incidence/Transmission", au["Incidence"]/au["Transmission"]];
    emit[tag, "Up.Reflection/Transmission", au["Reflection"]/au["Transmission"]];
    emit[tag, "Up.R(rprime)/Transmission", Rup[rprime]/au["Transmission"]];
    emit[tag, "Up.dR(rprime)/Transmission", Rup'[rprime]/au["Transmission"]];
    emit[tag, "seconds", Round[AbsoluteTime[] - t0, 0.1]];
  ], Print["ABORTED|", ToString[InputForm[c]]]],
  {c, cases}
];
