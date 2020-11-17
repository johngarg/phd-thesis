styles = Style[#, FontSize -> 18] & /@ {
    "\!\(\*SubscriptBox[\(\[Nu]\), \(1\)]\)", 
    "\!\(\*SubscriptBox[\(\[Nu]\), \(2\)]\)", 
    "\!\(\*SubscriptBox[\(\[Nu]\), \(3\)]\)", 
    "e", "\[Mu]", "\[Tau]", "d", "s", "b", "u", "c", "t"
};

fig = Show[
 ListLogLinearPlot[{
   {{0.0005, 1}, {0.005, 2}, {0.07, 3}}, (* neutrinos *)
   {{511 10^3, 1 - 0.1}, {106 10^6, 2 - 0.1}, {1.78 10^9, 3 - 0.1}}, (* charged leptons *)
   {{7 10^6, 1}, {120 10^6, 2}, {4.3 10^9, 3}}, (* down quarks *)
   {{3 10^6, 1 + 0.1}, {1.2 10^9, 2 + 0.1}, {174 10^9, 3 + 0.1}} (* up quarks *)
   },
  (*PlotStyle\[Rule]colors,*)
  PlotTheme -> "Detailed",
  PlotStyle -> {PointSize[0.03]},
  PlotMarkers -> {Style["\[Nu]?", FontSize -> 18], 
    Style["e", FontSize -> 18], Style["d", FontSize -> 18], 
    Style["u", FontSize -> 18]},
  AspectRatio -> 1/2,
  PlotRange -> {{10^-4, 10^12}, {0, 4}},
  GridLines -> {Automatic, {1, 2, 3}},
  FrameTicks -> {Automatic, {1, 2, 3}},
  ImageSize -> Large,
  FrameLabel -> {"Mass [eV]", "Generation"},
  LabelStyle -> Directive[Large, Black],
  Epilog -> {Rotate[
     Inset[Style["\!\(\*SubscriptBox[\(m\), \(\[Nu]\)]\) < 0.12 eV", 
       FontSize -> 18, 
       FontColor -> Gray], {-1, 
       3}], \[Pi]/2],
    Rotate[
     Inset[Style["\!\(\*SubscriptBox[\(m\), \(max\)]\) > 0.05 eV", 
       FontSize -> 18, 
       FontColor -> ColorData[97, "ColorList"][[7]]], {-3.8, 
       1.1}], \[Pi]/2]}
  ],
 ListLogLinearPlot[{{0.12, 0}, {0.12, 4}}, Joined -> True, 
  PlotStyle -> {Dashed, Gray}],
 ListLogLinearPlot[{{0.05, 0}, {0.05, 4}}, Joined -> True, 
  PlotStyle -> {Dashed, ColorData[97, "ColorList"][[7]]}]
 ];

Export["../../img/chapter_1/fermion_mass_plot.pdf", fig];


