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

def main(starting_time, offset_in_minutes)
  starting_time = 0 unless starting_time
  offset_in_minutes = 0 unless offset_in_minutes
  results_string = "results_start_" + starting_time.to_s + "_offset_" + offset_in_minutes.to_s
  results = Pathname.new(Dir.pwd) + "Data" + results_string
  puts results
  
  # Delete existing directory -- make sure nothing valuable is stored here
  FileUtils.remove_dir(results, force = true)
  FileUtils.mkdir(results)
  
  
  # path = Pathname.new(Dir.pwd) # gives current directory
  puts Dir.pwd
  path = Pathname.new(Dir.pwd) + "CWAMEnsembles"
  
  id = 0
  sum = 0
  test_sum = 0
  path.each_child do |folder|
    current_dir = path + folder
    puts current_dir
    children = current_dir.children
    time_upper_bound = parse_time_from_filename(children.first) + offset_in_minutes * SECONDS_PER_MINUTE + starting_time * SECONDS_PER_MINUTE
    time_lower_bound = parse_time_from_filename(children.first) + starting_time * SECONDS_PER_MINUTE
    children_in_offset = 0
    children.each do |child|
      time = parse_time_from_filename(child)
      break unless time <= time_upper_bound
      children_in_offset += 1 unless time < time_lower_bound
    end
    current_dir.each_child do |file_name|
      time = parse_time_from_filename(file_name)
      if time_lower_bound <= time and time <= time_upper_bound
        writeable_name = Pathname.new(results) + (id.to_s + ".dat")
        print time, "  ", writeable_name, "\n"
        w = File.new(writeable_name, 'a')
        temp = File.new(file_name, 'r')
        e = temp.each_line
        w.write(e.next)
        line = e.next.split
        w.write(line[0] + " ")
        n = line[1].to_f / children_in_offset
        test_sum += line[1].to_f
        sum += n
        #puts line[1].to_f, n
        w.write(n.to_s + "\n")
        begin
          loop do 
            w.write(e.next)
          end
        rescue StopIteration
          break
        ensure
          w.close
          temp.close
        end
        
        id += 1
      end
    end
  end
  puts "Total Probability (should be close to 1): " + sum.to_s
end

if __FILE__ == $0
  require './create_input'
  if ARGV.length % 2 == 1
    puts "Must provide one flag for each of the two possible parameters."
    puts "eg: ruby convert_files.rb -s 15 -o 30 -dthresh 0.8 -nethresh 0.8"
  else
    s = false
    o = false
    deviation_threshold = false
    node_edge_threshold = false
    ARGV.each do |arg|
      # The next two lines need ``== true`` because they are being used
      # as flags for themselves! Any non false / nil var will be true-like in an if statement.
      s = arg if s == true
      o = arg if o == true
      deviation_threshold = arg if deviation_threshold == true
      node_edge_threshold = arg if node_edge_threshold == true
      if arg == "-s"
        s = true
      elsif arg == "-o"
        o = true
      elsif arg == "-dthresh"
        deviation_threshold = true
      elsif arg == "-nethresh"
        node_edge_threshold = true
      end
    end
    main(s, o)
    create_input(deviation_threshold, node_edge_threshold)
  end
end