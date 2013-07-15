# Pass it params: single weather file or ensemble (of k)
require 'pathname'
require 'fileutils'

if __FILE__ == $0

  seed = 1
  number_of_members = 1
  num_weather_points = 1000
  min_lat = 38
  max_lat = 45.5
  min_lon = -93
  max_lon = -83
  max_prob = 1.0
  min_prob = 0.7
  altitudes = ["10000", "15000", "20000", "25000", "30000"]
  file_probabilities = { 1 => 0.8, 2 => 0.1, 3 => 0.02, 4 => 0.02, 5 => 0.02, 6 => 0.02 }
  weather_dir = "Data"

  ARGV.each_with_index do |arg, index|
    seed                = ARGV[index+1].to_i if arg == "-seed"
    number_of_members   = ARGV[index+1].to_i if arg == "-num_members"
    num_weather_points  = ARGV[index+1].to_i if arg == "-num_weather_points"
    max_lat             = ARGV[index+1].to_f if arg == "-max_lat"
    min_lat             = ARGV[index+1].to_f if arg == "-min_lat"
    max_prob            = ARGV[index+1].to_f if arg == "-max_prob"
    min_prob            = ARGV[index+1].to_f if arg == "-min_prob"
    weather_dir         = ARGV[index+1]      if arg == "-weather_dir"
  end

  begin
    FileUtils.mkdir(weather_dir)
  rescue Errno::EEXIST
  end

  seed_specific_folder = "seed" + seed.to_s + "_numMembers" + number_of_members.to_s
  seed_dir = Pathname.new(weather_dir) + seed_specific_folder
  begin
    FileUtils.mkdir(seed_dir)
  rescue Errno::EEXIST
  end

  1.upto(number_of_members) do |member_number|
    probability_written = 0
    file = "Seed_" + seed.to_s + "_20090618T060000_Member" + member_number.to_s + ".dat"
    rng = Random.new(seed)
    file_name = Pathname.new(seed_dir) + file

    File.open(file_name, 'w') do |f|
      f.write("Ensemble Member: ")
      f.write(member_number.to_s)
      f.write("/")
      f.write(number_of_members.to_s + "\n")
      f.write("Probability ")
      if member_number < number_of_members
        f.write(file_probabilities[member_number].to_s) 
        probability_written += file_probabilities[member_number]
      else if member_number == number_of_members
        f.write((1 - probability_written).to_s)
      end
      f.write("\n")
      
      num_weather_points.times do
        # Generate latitude between 38 and 45.5
        # Generate longitude between -93 and -83
        lat = rng.rand * (max_lat - min_lat) + min_lat
        lon = rng.rand * (max_lon - min_lon) + min_lon
        prob = rng.rand * (max_prob - min_prob) + min_prob

        altitudes.each do |alt|
          f.write("300,0,0,")
          f.write(("%.6f" % lat).to_s + "," + ("%.6f" % lon).to_s + "," + alt)
          f.write("," + ("%.3f" % prob).to_s + "\n")
        end
      end
    end
  end
  end


  require './create_input'
  if false 
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

    create_input(dshift, ddrop, angle, deviation_threshold, node_edge_threshold, output_name, seed_dir, c_input_file, weather_cell_width, quadrant_size, lane_width, max_fix_nodes, oper_flex)
  end
end
