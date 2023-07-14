(* ::Package:: *)

(* ::Section:: *)
(*Raw Data*)


NotebookDirectory[]


windfile = "/Users/gagebonner/Downloads/KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_25_DES_1688757095019.nc";
(*https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-WIND-PUM-012-002-005.pdf*)
waterTemperatureFile="/Users/gagebonner/Downloads/global-reanalysis-phy-001-031-grepv2-daily_1688755986955.nc";
(*https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-031.pdf*)
nutrientsFile="/Users/gagebonner/Downloads/cmems_mod_glo_bgc_my_0.25_P1D-m_1688755326762.nc";
(*https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-029.pdf*)


Import[windfile, "Dimensions"]
Import[waterTemperatureFile, "Dimensions"]
Import[nutrientsFile, "Dimensions"]


(* ::Section:: *)
(*Wind*)


(* ::Subsection:: *)
(*Longitude*)


(*Longitude is given in degrees East; subtract 360 to convert to degrees E/W.*)


windlon=Normal[Import[windfile,"/lon"]] - 360;
windlon[[1;;10]]


(* ::Subsection:: *)
(*Latitude*)


(*Longitude is given in degrees N/S*)


windlat=Normal[Import[windfile,"/lat"]];
windlat[[1;;10]]


(* ::Subsection:: *)
(*Time*)


(*Time is given in seconds since 1990-01-01 00:00:00*)
(*The function AbsoluteTime gives the number of seconds since January 1, 1900.*)
(*The function DayCount[DateObject[{1,1,1,0,0,0}], t] converts t to a Rata Die time (the number of days since January 1, 1 A.D.)*)


windtime=Normal[Import[windfile,"/time"]];
dateobjects=Map[t|->FromAbsoluteTime[AbsoluteTime[{1990,1,1,0,0,0}]+t],windtime];
windtime=Map[t|->DayCount[DateObject[{1,1,1,0,0,0}], t], dateobjects];
windtime[[1;;10]]


(* ::Subsection:: *)
(*Components*)


(*u = /eastward_wind, v = /northward_wind, convention is uv[time, lat, lon]*)
(*Have eastward and westward winds given in m/s; multiply by 86.4 to convert to km/day*)
(*Values must be scaled by 0.01*)
(*Missing/NaNs are -32767*)


WindProcess[x_]:=If[x==-32767,0,x*0.01*86.4]


uwind=Normal[Import[windfile,"/eastward_wind"]];
uwind=Map[WindProcess,uwind,{3}];
uwind=Transpose[uwind, Ordering[{3,2,1}]];

vwind=Normal[Import[windfile,"/northward_wind"]];
vwind=Map[WindProcess,vwind,{3}];
vwind=Transpose[vwind, Ordering[{3,2,1}]];


{Length[windlon],Length[windlat],Length[windtime]}
Dimensions[uwind]
Dimensions[vwind]


(* ::Subsection:: *)
(*Exporting*)


Export[NotebookDirectory[]<>"wind_glor.mat",
	{
	"lon"->N[windlon], 
	"lat"->N[windlat], 
	"time"->N[windtime], 
	"u"->N[uwind], 
	"v"->N[vwind]
	}
]


(* ::Section:: *)
(*Water*)


(* ::Subsection:: *)
(*Longitude*)


(*Longitude is given in degrees E/W*)


waterlon=Normal[Import[waterTemperatureFile,"/longitude"]];
waterlon[[1;;10]]


(* ::Subsection:: *)
(*Latitude*)


(*Longitude is given in degrees N/S*)


waterlat=Normal[Import[waterTemperatureFile,"/latitude"]];
waterlat[[1;;10]]


(* ::Subsection:: *)
(*Time*)


(*Time is given in days since 1950-01-01 00:00:00*)


watertime=Normal[Import[waterTemperatureFile,"/time"]];
watertime=Map[t|->DayCount[DateObject[{1,1,1,0,0,0}], DateObject[{1950,1,1,0,0,0}]]+t,watertime];
watertime[[1;;10]]


(* ::Subsection:: *)
(*Components*)


(*u = /uo_glor, v = /vo_glor, convention is uv[time, depth, lat, lon]*)
(*Have eastward and westward winds given in m/s; multiply by 86.4 to convert to km/day*)
(*no scale factor*)
(*Missing/NaNs are 9.96921e+36f*)


WaterProcess[x_]:=If[x>10^30,0,x*1*86.4]


uwater=Normal[Import[waterTemperatureFile,"/uo_glor"]];
uwater=uwater[[All,1,All,All]]; (*extract depth axis*)
uwater=Map[WaterProcess,uwater,{3}];
uwater=Transpose[uwater, Ordering[{3,2,1}]];

vwater=Normal[Import[waterTemperatureFile,"/vo_glor"]];
vwater=vwater[[All,1,All,All]]; (*extract depth axis*)
vwater=Map[WaterProcess,vwater,{3}];
vwater=Transpose[vwater, Ordering[{3,2,1}]];


{Length[waterlon],Length[waterlat],Length[watertime]}
Dimensions[uwater]
Dimensions[vwater]


(* ::Subsection:: *)
(*Exporting*)


Export[NotebookDirectory[]<>"water_glor.mat",
	{
	"lon"->N[waterlon], 
	"lat"->N[waterlat], 
	"time"->N[watertime], 
	"u"->N[uwater], 
	"v"->N[vwater]
	}
]


(* ::Section:: *)
(*Temperature*)


(* ::Subsection:: *)
(*Longitude*)


(*Longitude is given in degrees E/W*)


templon=Normal[Import[waterTemperatureFile,"/longitude"]];
templon[[1;;10]]


(* ::Subsection:: *)
(*Latitude*)


(*Longitude is given in degrees N/S*)


templat=Normal[Import[waterTemperatureFile,"/latitude"]];
templat[[1;;10]]


(* ::Subsection:: *)
(*Time*)


(*Time is given in days since 1950-01-01 00:00:00*)


temptime=Normal[Import[waterTemperatureFile,"/time"]];
temptime=Map[t|->DayCount[DateObject[{1,1,1,0,0,0}], DateObject[{1950,1,1,0,0,0}]]+t,temptime];
temptime[[1;;10]]


(* ::Subsection:: *)
(*Components*)


(*temp = /thetao_glor, convention is temp[time, depth, lat, lon]*)
(*Temperature is given in degrees Celcius*)
(*no scale factor*)
(*Missing/NaNs are 9.96921e+36f*)


TemperatureProcess[x_]:=If[x>10^30,0,x]


temp=Normal[Import[waterTemperatureFile,"/thetao_glor"]];
temp=temp[[All,1,All,All]]; (*extract depth axis*)
temp=Map[TemperatureProcess,temp,{3}];
temp=Transpose[temp, Ordering[{3,2,1}]];


{Length[templon],Length[templat],Length[temptime]}
Dimensions[temp]


(* ::Subsection:: *)
(*Exporting*)


Export[NotebookDirectory[]<>"temp_glor.mat",
	{
	"lon"->N[templon], 
	"lat"->N[templat], 
	"time"->N[temptime], 
	"u"->N[temp]
	}
]


(* ::Section:: *)
(*Nutrients*)


(* ::Subsection:: *)
(*Longitude*)


(*Longitude is given in degrees E/W*)


nutlon=Normal[Import[nutrientsFile,"/longitude"]];
nutlon[[1;;10]]


(* ::Subsection:: *)
(*Latitude*)


(*Longitude is given in degrees N/S*)


nutlat=Normal[Import[nutrientsFile,"/latitude"]];
nutlat[[1;;10]]


(* ::Subsection:: *)
(*Time*)


(*Time is given in hours since 1950-01-01 00:00:00*)


nuttime=Normal[Import[nutrientsFile,"/time"]];
nuttime=Map[t|->DayCount[DateObject[{1,1,1,0,0,0}], DateObject[{1950,1,1,0,0,0}]]+t/24,nuttime];
nuttime = nuttime - 1; (*Shifted by 1 wrt water/wind*)
nuttime[[1;;10]]


(* ::Subsection:: *)
(*Components*)


(*nitrogen = /no3, convention is no3[time, depth, lat, lon]*)
(*Nitrogen is given in mmol/m^3*)
(*no scale factor*)
(*Missing/NaNs are 9.96921e+36f*)


NitrogenProcess[x_]:=If[x>10^30,0,x]


nitrogen=Normal[Import[nutrientsFile,"/no3"]];
nitrogen=nitrogen[[All,1,All,All]]; (*extract depth axis*)
nitrogen=Map[NitrogenProcess,nitrogen,{3}];
nitrogen=Transpose[nitrogen, Ordering[{3,2,1}]];


{Length[nutlon],Length[nutlat],Length[nuttime]}
Dimensions[nitrogen]


(* ::Subsection:: *)
(*Exporting*)


Export[NotebookDirectory[]<>"no3_glor.mat",
	{
	"lon"->N[nutlon], 
	"lat"->N[nutlat], 
	"time"->N[nuttime], 
	"u"->N[nitrogen]
	}
]
