pmns = Block[{
        \[Theta]12 = 33.44 Degree, 
        \[Theta]23 = 49.0 Degree, 
        \[Theta]13 = 8.57 Degree, 
        \[Delta] = 195 Degree, 
        c12, c13, s13, s23, c23, s12
    },

    c12 = Cos[\[Theta]12];
    c13 = Cos[\[Theta]13];
    s13 = Sin[\[Theta]13];
    s23 = Sin[\[Theta]23];
    c23 = Cos[\[Theta]23];
    s12 = Sin[\[Theta]12];

    {
        {c12 c13,                                s12 c13,                                s13 Exp[-I \[Delta]]},
        {-s12 c23 - c12 s23 s13 Exp[I \[Delta]], c12 c23 - s12 s23 s13 Exp[I \[Delta]],  s23 c13},
        {s12 c23 - c12 c23 s13 Exp[I \[Delta]],  -c12 s23 - s12 c23 s13 Exp[I \[Delta]], c23 c13}
    }
 ];

 absPmns = Map[Abs, pmns, -1];
 Print["Magnitude of PMNS matrix elements:"];
 unit = "em"
 Table[
     Table[
         Print["\\rule{"<> ToString[i[[j]]] <> unit <> "}{" <> ToString[i[[j]]] <> unit <> "}" <> If[j === 3, " //", " &"]],
         {j, Length[i]}
     ]; 
     Print[], 
 {i, absPmns}
 ];
