(* ::Package:: *)

(* ::Subsection:: *)
(*Units*)


(* ::Subsubsection:: *)
(*Messages*)


MyUnits::unknown="Unknown units.  Cannot convert.";


MyUnits::incomp="Incompatible units.  Supplemented natural units.";


(* ::Subsubsection:: *)
(*Conversion*)


UnitsTypes={"Energy","Time","Length","Area","Volume","Mass","Temperature","Power","Force"};


UnitsNames={
{"eV","keV","MeV","GeV","TeV","PeV","joule","erg","calorie"},
{"fsec","psec","nsec","\[Mu]sec","msec","sec","ksec","Msec","Gsec","Tsec","minute","hour","day","week","year"},
{"angstrom","fm","pm","nm","\[Mu]m","mm","cm","meter","km","AU","lightyear","pc","kpc","Mpc","Gpc","mile","foot","inch","yard"},
{"\!\(\*SuperscriptBox[\"cm\", \"2\"]\)","\!\(\*SuperscriptBox[\"meter\", \"2\"]\)","\!\(\*SuperscriptBox[\"km\", \"2\"]\)","acre","hectare","ab","fb","pb","nb","\[Mu]b","mb","barn"},
{"\!\(\*SuperscriptBox[\"cm\", \"3\"]\)","\!\(\*SuperscriptBox[\"meter\", \"3\"]\)","\!\(\*SuperscriptBox[\"km\", \"3\"]\)","quart","liter","pint","cup","gallon"},
{"AMU","gram","kg","ton","pound","ounce"},
{"kelvin","celsius","fahrenheit"},
{"GeV \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"1\"}]]\)","watt","kg \!\(\*SuperscriptBox[\"meter\", \"2\"]\) \!\(\*SuperscriptBox[\"sec\", \"3\"]\)"},
{"newton","kg cm \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"2\"}]]\)","kg meter \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"2\"}]]\)","dyne","pound","kpound"}
};


UnitsValueCKS={
{10^-9,10^-6,10^-3,1,10^3,10^6,6.24150974*10^9,624.150974,2.61144768*10^10},
{10^-15,10^-12,10^-9,10^-6,10^-3,1,10^3,10^6,10^9,10^12,60,60 60,24 60 60,7 24 60 60,24 60 60 365.242199},
{10^-8,10^-13,10^-10,10^-7,10^-4,10^-1,1,10^2,10^5,1.49598*10^13,9.4605284*10^17,3.08568025*10^18,3.08568025*10^21,3.08568025*10^24,3.08568025*10^27,160934.4,30.48,2.54,91.44},
{1,10^4,10^10,40468564.2,10^8,10^-42,10^-39,10^-36,10^-33,10^-30,10^-27,10^-24},
{1,10^6,10^15,946.352946,1000,473.176473,236.588236,3.7854118 10^3},
{1.660538782*10^-27,10^-3,1,10^3,0.45359237,0.0283495231},
{1,274.15,255.927778},
{1,6.24150974*10^9,6.24150974*10^9},
{100,1,100,10^-3,444.822162,444.822162 10^3}
};


UnitsNaturalCKS={GeV,1.51926778 * 10^24 GeV^-1,5.06773182*10^13 GeV^-1,(5.06773182*10^13 GeV^-1)^2,(5.06773182*10^13 GeV^-1)^3,5.60958921*10^26 GeV,8.6173423*10^-14 GeV,6.58211814*10^-25 GeV^2,1.231618*10^-8 GeV^2};


UnitsStdList={
{"CKS",{"GeV","sec","cm","\!\(\*SuperscriptBox[\"cm\", \"2\"]\)","\!\(\*SuperscriptBox[\"cm\", \"3\"]\)","kg","kelvin","GeV \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"1\"}]]\)","kg cm \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"2\"}]]\)"}},
{"MKS",{"GeV","sec","meter","\!\(\*SuperscriptBox[\"meter\", \"2\"]\)","\!\(\*SuperscriptBox[\"meter\", \"3\"]\)","kg","kelvin","GeV \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"1\"}]]\)","kg meter \!\(\*SuperscriptBox[\"sec\", 
RowBox[{\"-\", \"2\"}]]\)"}}
};


UnitsStd=UnitsStdList[[1,2]];
UnitsNatural=UnitsNaturalCKS;
UnitsValue=UnitsValueCKS;


(* ::Subsubsection:: *)
(*Constants*)


ConstTypes={"General","SM","Atomic","Astrophysics","Cosmology"};


ConstValues={
{{"c",29979245800 cm sec^-1},
{"hbar",6.58211814*10^-25 GeV sec},
{"kB",8.6173423*10^-14 GeV kelvin^-1},
{"Mpl",(6.70881 10^-39)^(-1/2) GeV}},
{{"melectron",0.510998910 MeV},
{"mmuon",105.658369 MeV},
{"mtau",1776.82 MeV},
{"mup",2.49 MeV},
{"mdown", 5.05 MeV},
{"mstrange",101 MeV},
{"mcharm",1.27 GeV},
{"mbottom",4.19 GeV},
{"mtop", 172 GeV},    
{"mproton",938.272013 MeV},
{"mneutron",939.565346 MeV},
{"mpipm", 139.570 MeV},
{"mpizero",134.977 MeV},
{"mW",80.398 GeV},
{"mZ",91.1876 GeV},
{"GF",1.16637 10^-5 GeV^-2},
{"\[Alpha]EM",7.2973525698 10^-3},
{"\[Alpha]1",7.2973525698 10^-3  /(1-0.23116)},
{"\[Alpha]2",7.2973525698 10^-3 / 0.23116},
{"\[Alpha]3",0.1184 },
{"sW2",0.23116},
{"vHiggs",246 GeV}},
{{"NA",6.02214179 10^23  mol^-1},
{"relectron",2.817940325 fm}},
{{"Msun", 1.9884 10^30  kg},
{"Rsun",6.9551 10^8 meter},
{"Lsun",3.8427 watt},
{"Mearth",5.9722 10^24 kg},
{"Rearth",6.378137 10^6 meter}},
{{"\[Rho]c",1.05368 10^-5 h^2 GeV cm^-3},
{"s0",2889.2 cm^-3},
{"H0",100 h km sec^-1 Mpc^-1},
{"h0",0.673},
{"\[CapitalOmega]cdm0",0.12038 h^-2},
{"\[CapitalOmega]m0",0.3183},
{"\[CapitalOmega]b0",0.022032 h^-2},
{"\[CapitalOmega]\[Gamma]0",2.471 10^-5 h^-2},
{"\[CapitalOmega]\[CapitalLambda]0",0.6817},
{"\[Eta]baryon",6.23 10^-10},
{"n\[Gamma]0",410.5 cm^-3},
{"T0",2.7255 kelvin},
{"zeq",3391}}
};


(* ::Subsubsection:: *)
(*Functions*)


PrintUnits[]:=Grid[Table[{Grid[{{UnitsTypes[[i]]},{Style[StringJoin["[",UnitsStd[[i]],"]"],Red]}}],Grid[Transpose[{UnitsNames[[i]],UnitsValue[[i]]}//N],Background->{None,{{LightYellow,LightGray}}},Alignment->{Left},Spacings->{3},ItemSize->9]},{i,1,Length[UnitsTypes]}],Frame->All,Background->{LightBlue},Alignment->{Center}]


ClearUnits[]:=Module[{},
Apply[Clear,Flatten[UnitsNames]];
$Assumptions=Map[ToExpression[#]>0&,Flatten[UnitsNames]];
]


RemoveUnits[exp_]:=exp/.Map[ToExpression[#]->1&,Flatten[UnitsNames]]


NaturalRemoveUnits[exp_]:=RemoveUnits[NaturalUnits[exp]]


SetStdUnits[name_]:=Module[{pos},

If[Length[Position[UnitsStdList,name]]==0,Message[MyUnits::unknown];Return[];];

pos=Position[UnitsStdList,name][[1,1]];
UnitsStd=UnitsStdList[[pos,2]];

UnitsValue=Table[UnitsValueCKS[[i]]/UnitsValueCKS[[i,Position[UnitsNames,UnitsStd[[i]]][[1,2]]]],{i,1,Length[UnitsStd]}]//N;

UnitsNatural=Table[UnitsNaturalCKS[[i]]*UnitsValueCKS[[i,Position[UnitsNames,UnitsStd[[i]]][[1,2]]]],{i,1,Length[UnitsStd]}]//N;


]


StdUnits[exp_]:=Module[{replace},
replace=Flatten[Table[ToExpression[UnitsNames[[i,j]]]->ToExpression[UnitsStd[[i]]]UnitsValue[[i,j]], {i,1,Length[UnitsNames]},{j,1,Length[UnitsNames[[i]]]}]];


exp/.Select[replace,Length[Variables[#[[1]]]]==1&]//N//Simplify
]


NaturalUnits[exp_]:=Module[{replace},

replace=Table[ToExpression[UnitsStd[[i]]]->UnitsNatural[[i]], {i,1,Length[UnitsStd]}];

StdUnits[exp]/.Select[replace,Length[Variables[#[[1]]]]==1&]//Simplify

]


InUnits[exp_,units_]:=Module[{unitsvars,expvars},
unitsvars=Variables[units];
expvars=Variables[NaturalUnits[exp/units]];


If[Length[Select[Map[Position[UnitsNames,ToString[#]]&,unitsvars],Length[#]==0&]]!=0,Message[MyUnits::unknown];Return[exp];];
If[Length[expvars]!=0,Message[MyUnits::incomp];];

NaturalUnits[exp/units] units //N
]


InUnits[units_]:=InUnits[#,units]&


ClearConsts[]:=Module[{vars},
vars=Flatten[ConstValues,1][[1;;,1]];

Apply[Clear,vars];
]


SetConsts[]:=Module[{},
ClearConsts[];
Map[Apply[Set,{ToExpression[#[[1]]],#[[2]]}]&,Flatten[ConstValues,1]];
]


PrintConsts[]:=Grid[Table[{ConstTypes[[i]],Grid[ConstValues[[i]]//N,Background->{None,{{LightYellow,LightGray}}},Alignment->{Left},Spacings->{3},ItemSize->10]},{i,1,Length[ConstTypes]}],Frame->All,Background->{LightBlue},Alignment->{Center}]
