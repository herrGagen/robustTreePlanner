@echo off

mkdir XML
setlocal enabledelayedexpansion

:: n is node-edge threshold varying from .2 to .9
:: d is weather day 	For now, we are using starting time as a standin for different days (do we have different days?)
:: o is weather offset, we're testing 0, 30, 60, 90, 120 minutes here.

for /l %%n in (9, -1, 2) do @(
	for %%d in (6 7 8) do @(
		for /l %%o in (0,30,120) do @(
			set weatherDay=201906180600/0%%d00
			set outputName=XML/neThresh_0p%%n_weatherDay_%%d_offset_%%o
			call run.bat -nethresh 0.%%n -iname !weatherDay! -o %%o -oname !outputName! -cellwidth .1 -lanewidth 3 -angle 270 -s 0 -twname tempWeather -cinput originalInput.txt -quadrantsize 90 -operflex 1.0 2.0 4.0 6 8 10 12 14 16 18 20 22 24 -fixnodes 1 -demandshift 1 -demanddrop 4 -dthresh 0.8 -concat
		)
	)
)

:: call compile_and_run.bat -twname basic_example -cinput originalInput.txt -iname 0600 -quadrantsize 90 -operflex 1.0 2.0 4.0 8.0 -fixnodes 1 -demandshift 1 -demanddrop 4 -s 0 -o 45 -angle 235 -dthresh 0.8 -nethresh 0.7 -cellwidth 2 -lanewidth 3 -oname q4d80n70cw2lw3
:: call compile_and_run.bat -twname basic_example -cinput originalInput.txt -iname 0600 -quadrantsize 90 -operflex 1.0 2.0 4.0 8.0 -fixnodes 1 -demandshift 1 -demanddrop 6 -s 0 -o 15 -angle 235 -dthresh 0.8 -nethresh 0.7 -cellwidth 2 -lanewidth 3 -oname q4d80n70cw2lw3 -concat