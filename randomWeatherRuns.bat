@echo off
setlocal enabledelayedexpansion
set NUM_MEMBERS=1
set MIN_NUMEL=100
set MAX_NUMEL=1000
set MIN_SEED=1
set MAX_SEED=1000
set xmlName=blank

FOR /L %%n IN (%MIN_NUMEL%,100,%MAX_NUMEL%) DO @(
  FOR /L %%s IN (1,1,5) DO @(
    FOR /L %%a IN (0, 90, 270) DO @(
      FOR %%w IN (-1 2 4) DO @(
        set xmlName=seed%%s_angle%%a_laneWidth%%w_numPoints%%n_numMembers%NUM_MEMBERS%
        ruby generate_random_weather.rb -seed %%s -lanewidth %%w -angle %%a -operflex .82 2 4 -num_weather_points %%n -num_members %NUM_MEMBERS% -oname !xmlName!
        RobustTree.exe -cinput inputs.txt
      )
    )
  )
)