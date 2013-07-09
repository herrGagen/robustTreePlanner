# Takes multiple files and merges them into one directory.  
# It will also adjust the probabilities of each time to accomodate multiple times

# BEFORE RUNNING MAKE SURE RESULTS IS BACKED UP (It will get deleted)

require 'fileutils'
require 'pathname'
require 'date'

SECONDS_PER_MINUTE = 60

# takes in a filename, strips out all non-numerals, 
# and passes it to DateTime to parse.  
# Returns an instance of the Time class.
def parse_time_from_filename(filename)
  DateTime.parse(filename.basename.to_s.gsub(/[^0-9]/, "")).to_time
end

def main(starting_time, offset_in_minutes, results_string=false, weather_dir=false)
  starting_time       = 0 unless starting_time
  offset_in_minutes   = 0 unless offset_in_minutes
  results_string      = "results_start_" + starting_time.to_s + "_offset_" + offset_in_minutes.to_s unless results_string
  results             = Pathname.new(Dir.pwd) + "Data" + results_string
  puts results
  
  weather_dir         = "CWAMEnsembles" unless weather_dir
  
  offset_in_minutes = offset_in_minutes.to_i
  starting_time = starting_time.to_i
  
  # Delete existing directory -- make sure nothing valuable is stored here
  FileUtils.remove_dir(results, force = true)
  FileUtils.mkdir(results)
  
  
  # path = Pathname.new(Dir.pwd) # gives current directory
  puts Dir.pwd
  path = Pathname.new(Dir.pwd) + weather_dir
  
  min_time = Time.new(Float::MAX)
  max_time = Time.new(Float::MIN)

  # First pass determines the range of times files can have
  # Could make this more efficient by only checking files from one Ensemble (all Member1's, for instance)
  path.each_child do |file_name|
    begin
      name_ary = file_name.basename.to_s.split("_")
      time = DateTime.parse(name_ary[2]).to_time # The 3rd part of the string contains the relevant time
      min_time = time if time < min_time
      max_time = time if time > max_time
    rescue TypeError
      print "Error parsing a file: ", file_name, "\n"
    end
  end
  
  raise "Files not in expected location or the filenames do not conform to ``DevProb_time1_time2_MemberX.dat``" if min_time == Time.new(Float::MAX) or max_time == Time.new(Float::MIN)
  time_upper_bound = min_time + offset_in_minutes * SECONDS_PER_MINUTE + starting_time * SECONDS_PER_MINUTE
  time_lower_bound = min_time + starting_time * SECONDS_PER_MINUTE

  
  # Second pass determines how many members of each Ensemble are within the time constraints.
  children_in_offset = 0
  member_name = path.children.first.basename.to_s.split("_")[3]
  path.each_child do |file_name|
    name_ary = file_name.basename.to_s.split("_")
    if name_ary[3] == member_name
      time = DateTime.parse(name_ary[2]).to_time
      if time <= time_upper_bound and time >= time_lower_bound
        children_in_offset += 1
      end
    end
  end

  print "Children in offset: ", children_in_offset, "\n" # I'm worried this value is incorrect  

  # Third pass writes files to a format readable by the C code
  sum = 0
  id = 0
  path.each_child do |file_name|
    name_ary = file_name.basename.to_s.split("_")
#     print file_name
    time = DateTime.parse(name_ary[2]).to_time # The 3rd part of the string contains the relevant time
    if time <= time_upper_bound and time >= time_lower_bound
      writeable_name = Pathname.new(results) + (id.to_s + ".dat")
#       print time, "  ", writeable_name, "\n"
      w = File.new(writeable_name, 'a')
      temp = File.new(file_name, 'r')
      begin
        e = temp.each_line
      
        # Need to find the Ensemble member number and  probability to write out
        4.times { e.next } # Burn through first 4 rows of the ensemble file (we don't need them)
        line = e.next
        line.slice!(0, 2)
        # print "LINE: ", line, "\n"
        w.write("Ensemble " + line)
        line = e.next.split.last.to_f / children_in_offset
        sum += line
      
        # Will changed 2013-05-06
        # Made this output with a force to output 9 decimals, since before it would use exponential notation.
        w.write("Probability " + ("%.9f" % line).to_s + "\r\n")
        # w.write("Probability " + line.to_s + "\r\n") 

        line = e.next
        while line.split.first == "#"
          line = e.next
        end
      rescue StopIteration
        print "Current File is not weather: ", file_name, "\n"
        break
      ensure
        w.close unless w.closed?
        File.delete(writeable_name)
        temp.close unless temp.closed?
      end 
      unless w.closed? or temp.closed?
        begin
          loop do
            w.write(line) if line.split(",").last.to_f >= 0.70
            # print line.split(","), "\n"
            line = e.next
          end
        rescue StopIteration
          break
        ensure
          w.close unless w.closed?
          temp.close unless temp.closed?
        end
        
        id += 1
      end
    end
  end  

  puts "Total Probability (should be close to 1): " + sum.to_s
end

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
      quadrant_size       = arg if quadrant_size        == true # this is called angle_offset in "create_input.rb"
      lane_width          = arg if lane_width	          == true
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
    print "Quadrant Size (angle offset):  ", quadrant_size,           "\n" if quadrant_size
    print "Deviation threshold:           ", deviation_threshold,     "\n" if deviation_threshold
    print "Node edge threshold:           ", node_edge_threshold,     "\n" if node_edge_threshold    
    print "Output file name:              ", output_name,             "\n" if output_name    
    print "Temp weather dir:              ", temp_weather_name,       "\n" if temp_weather_name    
    print "Input weather dir:             ", weather_dir,             "\n" if weather_dir
    print "Weather cell width:            ", weather_cell_width,      "\n" if weather_cell_width
    print "Lane width:                    ", lane_width,              "\n" if lane_width
    print "Max Number of Fix Nodes:       ", max_fix_nodes,           "\n" if max_fix_nodes
    print "Operational Flexibility:       ", oper_flex,               "\n" if oper_flex != []
    
    main(s, o, temp_weather_name, weather_dir)
    create_input(dshift, ddrop, angle, deviation_threshold, node_edge_threshold, output_name, temp_weather_name, c_input_file, weather_cell_width, quadrant_size, lane_width, max_fix_nodes, oper_flex)
  end
end
