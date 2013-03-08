require 'fileutils'
require 'pathname'
require 'date'

def create_input(dthresh, nethresh, oname)
  deviation_threshold = dthresh ? dthresh : 0.8
  node_edge_threshold = nethresh ? nethresh : 0.8
  output_name = oname ? oname : Time.now.to_s.slice(0, 19).gsub(" ", "_").gsub(":", "-")
  
  current_dir = Pathname.new(Dir.pwd)
  data_path = current_dir + "Data"
  weather_path = nil
  dp_children = data_path.children
  creation_times = dp_children.map { |e| File.ctime(e).to_i }
  weather_path = data_path.children[creation_times.index(creation_times.max)]
  
  input_file_path = current_dir + "inputs.txt"
  input = File.new(input_file_path, 'w')
  
  input.write("Data/demand.nom\n")
  input.write(weather_path.children.length.to_s + "\n")
  weather_path.each_child do |child|
    input.write((child.dirname.dirname.basename + child.dirname.basename + child.basename).to_s + "\n")
  end
  input.write(deviation_threshold.to_s + "\n")
  input.write(node_edge_threshold.to_s + "\n")
  input.write(output_name.to_s + "\n")
end
