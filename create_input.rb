require 'fileutils'
require 'pathname'
require 'date'

def create_input(a, dthresh, nethresh, oname, temp_weather_name, c_input_file)
  angle = a ? a : 0.0
  deviation_threshold = dthresh ? dthresh : 0.8
  node_edge_threshold = nethresh ? nethresh : 0.8
  output_name = oname ? oname : Time.now.to_s.slice(0, 19).gsub(" ", "_").gsub(":", "-")
  c_input_file = c_input_file ? c_input_file : "inputs.txt"
  
  current_dir = Pathname.new(Dir.pwd)
  data_path = current_dir + "Data"
  weather_path = nil
  dp_children = data_path.children
  creation_times = dp_children.map { |e| File.ctime(e).to_i }
  weather_path = temp_weather_name ? data_path + temp_weather_name : data_path.children[creation_times.index(creation_times.max)]
  print "Temp weather dir: ", weather_path # This may be redundant, but it indicates which one it selected
                                           # If there is no option set, it should use the most recently created folder
  
  input_file_path = current_dir + c_input_file
  input = File.new(input_file_path, 'w')
  
  input.write("Data/demand.nom\n")
  input.write(weather_path.children.length.to_s + "\n")
  weather_path.each_child do |child|
    input.write((child.dirname.dirname.basename + child.dirname.basename + child.basename).to_s + "\n")
  end
  input.write(angle.to_s + "\n")
  input.write(deviation_threshold.to_s + "\n")
  input.write(node_edge_threshold.to_s + "\n")
  input.write(output_name.to_s + "\n")
end
