(* ::Package:: *)

(* ::Section:: *)
(*Raw Data*)


windfile = "/Users/gagebonner/Downloads/KNMI-GLO-WIND_L3-REP-OBS_METOP-A_ASCAT_25_DES_1688757095019.nc";
waterTemperatureFile="/Users/gagebonner/Downloads/global-reanalysis-phy-001-031-grepv2-daily_1688755986955.nc";
nutrientsFile="/Users/gagebonner/Downloads/cmems_mod_glo_bgc_my_0.25_P1D-m_1688755326762.nc";


Import[windfile, "Dimensions"]
Import[waterTemperatureFile, "Dimensions"]
Import[nutrientsFile, "Dimensions"]


(* ::Section:: *)
(*Wind*)


(* ::Subsection:: *)
(*Longitude*)


(*Longitude is given in degrees; subtract 360 to convert to degrees E/W.*)


windlon=Normal[Import[windfile,"/lon"]] - 360;
windlon[[1;;10]]


(* ::Subsection:: *)
(*Latitude*)


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


(* ::Section:: *)
(*Water*)


(* ::Subsection:: *)
(*Longitude*)


waterlon=Normal[Import[waterTemperatureFile,"/longitude"]];
waterlon[[1;;10]]


(* ::Subsection:: *)
(*Latitude*)


waterlat=Normal[Import[waterTemperatureFile,"/latitude"]];
waterlat[[1;;10]]


(* ::Subsection:: *)
(*Time*)


(*Time is given in days since 1950-01-01 00:00:00*)


watertime=Normal[Import[waterTemperatureFile,"/time"]];
watertime=Map[t|->DayCount[DateObject[{1,1,1,0,0,0}], DateObject[{1950,1,1,0,0,0}]]+t,watertime];
watertime[[1;;10]]


(* ::Section:: *)
(*Nutrients*)


(* ::Subsection:: *)
(*Longitude*)


nutlon=Normal[Import[nutrientsFile,"/longitude"]];
nutlon[[1;;10]]


DatePlus[DateObject[{1950,1,1,0,0,0}],24806]


24806/365//N


DateObject[24806,CalendarType->"Gregorian"]


DateObject[CalendarType->"Julian"]


FromUnixTime[880934400]


DayCount[DateObject[{1,1,1,0,0,0}],DatePlus[DateObject[{1950,1,1,0,0,0}],24806]]
