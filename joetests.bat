@echo off

if not exist XML mkdir XML
setlocal enabledelayedexpansion

:: n is node-edge threshold varying from .2 to .9
:: d is weather day 	For now, we are using starting time as a standin for different days (do we have different days?)

echo Starting
:: o is weather offset, we're testing 0, 30, 60, 90, 120 minutes here.

:: First, only vary node edge threshold and day
:: edge threshold = .1 .2 .4 .8 .9
:: offset = 120
:: date = all days
:: start time = first starting time per weather set
for %%n in (1 2 8 9 10) do @(
	for /l %%o in (120,30,120) do @(
		for %%d in (20090424 20090427 20090618 20110722) do @(
      echo Third %%d
			set dirname=%cd%\Weather\%%d\*
			set firstTime=1
			for /f %%t in ('dir /d /b !dirname!') do @(
      echo Fourth %%t
				if !firstTime!==1 @(
					set firstTime=0
					echo set weatherDay=Weather\%%d\%%t
					echo !weatherDay!
					for /l %%a in (0,90,180,270) do @(
						set outputName=XML/neThresh_0p%%n_weatherTime_%%t_weatherDay_%%d_offset_%%o_angle_%%a
						if not exist !outputName!.xml @(							
							start /b cmd /c call run.bat -angle %%a -nethresh 0.%%n -iname !weatherDay! -o %%o -oname !outputName! -cellwidth 1 -lanewidth 3 -s 0 -twname tempWeather -cinput !outputName!.joetests.txt -quadrantsize 90 -operflex 1.0 2.0 4.0 6 8 10 12 14 16 18 20 22 24 -fixnodes 1 -demandshift 1 -demanddrop 4 -dthresh 0.8 -concat
						)
					)
				)
			)
		)
	)
)

:: Next, only vary duration and day
:: edge threshold = .9
:: offset = 0 30 60 90 120
:: date = all days
:: start time = first starting time per weather set
for %%n in (9) do @(
	for /l %%o in (0,30,120) do @(
		for %%d in (20090424 20090427 20090618 20110722) do @(
			set dirname=%cd%\Weather\%%d\*
			set firstTime=1
			for /f %%t in ('dir /d /b !dirname!') do @(
				if !firstTime!==1 @(
					set firstTime=0
					set weatherDay=Weather\%%d\%%t
					echo !weatherDay!
					for /l %%a in (0,90,180,270) do @(
						set outputName=XML/neThresh_0p%%n_weatherTime_%%t_weatherDay_%%d_offset_%%o_angle_%%a
						if not exist !outputName!.xml @(
							start /b cmd /c call run.bat -angle %%a -nethresh 0.%%n -iname !weatherDay! -o %%o -oname !outputName! -cellwidth 1 -lanewidth 3 -s 0 -twname tempWeather -cinput !outputName!.joetests.txt -quadrantsize 90 -operflex 1.0 2.0 4.0 6 8 10 12 14 16 18 20 22 24 -fixnodes 1 -demandshift 1 -demanddrop 4 -dthresh 0.8 -concat
						)
					)
				)			
			)
		)
	)
)

:: Finally loop over all weather combinations
:: edge threshold = .9 
:: offset = 120
:: date = all days
:: start time = all starting times
for %%n in (9) do @(
	for /l %%o in (120,30,120) do @(
		for %%d in (20090424 20090427 20090618 20110722) do @(
			set dirname=%cd%\Weather\%%d\*
			for /f %%t in ('dir /d /b !dirname!') do @(
				set weatherDay=Weather\%%d\%%t
				echo !weatherDay!
				for /l %%a in (0,90,180,270) do @(
					set outputName=XML/neThresh_0p%%n_weatherTime_%%t_weatherDay_%%d_offset_%%o_angle_%%a
					if not exist !outputName!.xml @(
						start /b cmd /c call run.bat -angle %%a -nethresh 0.%%n -iname !weatherDay! -o %%o -oname !outputName! -cellwidth 1 -lanewidth 3 -s 0 -twname tempWeather -cinput !outputName!.joetests.txt -quadrantsize 90 -operflex 1.0 2.0 4.0 6 8 10 12 14 16 18 20 22 24 -fixnodes 1 -demandshift 1 -demanddrop 4 -dthresh 0.8 -concat
					)
				)
			)		
		)
	)
)
