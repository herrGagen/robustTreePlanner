# Takes multiple files and merges them into one directory.  
# It will also adjust the probabilities of each time to accomodate multiple times

# BEFORE RUNNING MAKE SURE RESULTS IS BACKED UP

require 'fileutils'
require 'pathname'
require 'date'

SECONDS_PER_MINUTE = 60
offset_in_minutes = 15
starting_time = 15

def parse_time_from_filename(filename)
  # specific method that takes in a filename, strips out all non-numerals, 
  # and passes it to datetime to parse.  
  # Returns an instance of the Time class.
  DateTime.parse(filename.basename.to_s.gsub(/[^0-9]/, "")).to_time
end

results_string = "results_start_" + starting_time.to_s + "_offset_" + offset_in_minutes.to_s
results = Pathname.new(Dir.pwd) + "Data" + results_string
puts results
FileUtils.remove_dir(results, force = true)
# exit(0)
FileUtils.mkdir(results)


# path = Pathname.new(Dir.pwd) # gives current directory
puts Dir.pwd
path = Pathname.new(Dir.pwd) + "CWAMEnsembles"
# use dirname and basename, ie path.dirname, path.basename
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

puts "Total Probability: " + sum.to_s
