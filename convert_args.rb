# Takes multiple files and merges them into one directory.
# It will also adjust the probabilities of each time to accomodate multiple times

# BEFORE RUNNING MAKE SURE RESULTS IS BACKED UP (It will get deleted)

require 'fileutils'
require 'pathname'
require 'date'

if __FILE__ == $0
  require './create_input'
  if false # ARGV.length % 2 == 1
    puts "Must provide one flag for each of the two possible parameters."
    puts "eg: ruby convert_files.rb -s 15 -o 30 -dthresh 0.8 -nethresh 0.8"
  else
    s                   = false
    o                   = false
    dshift              = false
    ddrop               = false
    angle               = false
    deviation_threshold = false
    node_edge_threshold = false
    output_name         = false
    temp_weather_name   = false
    weather_dir         = false
    c_input_file        = false
    weather_cell_width  = false
    quadrant_size       = false
    lane_width 	        = false
    max_fix_nodes       = false

    ARGV.each do |arg|
      # The following lines need ``== true`` because they are being used
      # as flags for themselves! Any non false / nil var will be true-like in an if statement.
      s = arg if s == true
      o = arg if o == true
      dshift = arg if dshift == true #just boolean 0 or 1, whether or not we want to use demand shifting
      ddrop = arg if ddrop == true # same as demand shifting, 0 or 1 to indicate whether or not we do it
      angle               = arg if angle                == true
      deviation_threshold = arg if deviation_threshold  == true
      node_edge_threshold = arg if node_edge_threshold  == true
      output_name         = arg if output_name          == true
      temp_weather_name   = arg if temp_weather_name    == true
      weather_dir         = arg if weather_dir          == true
      c_input_file        = arg if c_input_file         == true
      weather_cell_width  = arg if weather_cell_width   == true
      quadrant_size       = arg if quadrant_size        == true # this is called angular_width in "create_input.rb"
      lane_width          = arg if lane_width	        == true
      max_fix_nodes       = arg if max_fix_nodes        == true # called num_fix_nodes

      if arg == "-s"
        s = true
      elsif arg == "-o"
        o = true
      elsif arg == "-angle"
        angle = true
      elsif arg == "-demandshift"
        dshift = true
      elsif arg == "-demanddrop"
        ddrop = true
      elsif arg == "-dthresh"
        deviation_threshold = true
      elsif arg == "-nethresh"
        node_edge_threshold = true
      elsif arg == "-oname"
        output_name = true
      elsif arg == "-twname"
        temp_weather_name = true
      elsif arg == "-iname"
        weather_dir = true
      elsif arg == "-cinput"
        c_input_file = true
      elsif arg == "-cellwidth"
        weather_cell_width = true
      elsif arg == "-quadrantsize"
        quadrant_size = true
      elsif arg == "-lanewidth"
        lane_width = true
      elsif arg == "-fixnodes"
        max_fix_nodes = true
      end
    end

    oper_flex = []
    oper_flex_flag = false
    ARGV.each do |arg|
      if arg == "-operflex"
        oper_flex_flag = true
      elsif arg.split("").first == "-"
        oper_flex_flag = false
      end
      if oper_flex_flag and arg.split("").first != "-"
        oper_flex << arg
      end
    end

    angle = (angle.to_f * Math::PI / 180).to_s if angle
    quadrant_size = (quadrant_size.to_f * Math::PI / 180).to_s if quadrant_size
    print "Angle (converted to radians):  ", angle,                   "\n" if angle
    print "Quadrant Size (angular width): ", quadrant_size,           "\n" if quadrant_size
    print "Deviation threshold:           ", deviation_threshold,     "\n" if deviation_threshold
    print "Node edge threshold:           ", node_edge_threshold,     "\n" if node_edge_threshold    
    print "Output file name:              ", output_name,             "\n" if output_name    
    print "Temp weather dir:              ", temp_weather_name,       "\n" if temp_weather_name    
    print "Input weather dir:             ", weather_dir,             "\n" if weather_dir
    print "Weather cell width:            ", weather_cell_width,      "\n" if weather_cell_width
    print "Lane width:                    ", lane_width,              "\n" if lane_width
    print "Max Number of Fix Nodes:       ", max_fix_nodes,           "\n" if max_fix_nodes
    print "Operational Flexibility:       ", oper_flex,               "\n" if oper_flex != []

    create_input(dshift, ddrop, angle, deviation_threshold, node_edge_threshold, output_name, temp_weather_name, c_input_file, weather_cell_width, quadrant_size, lane_width, max_fix_nodes, oper_flex)
  end
end
