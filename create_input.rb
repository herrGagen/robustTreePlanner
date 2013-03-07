require 'fileutils'
require 'pathname'


current_dir = Pathname.new(Dir.pwd)
data_path = current_dir + "Data"
for child in data_path.each_child
  # Use a for loop to make weather_path a var outside of the loop
  # Need to clean this up for general times
  weather_path = child if child.basename.to_s[0] == "r"
end

input_file_path = current_dir + "input.txt"
input = File.new(input_file_path, 'a')

input.write("Data/demand.nom\n")
input.write(weather_path.children.length.to_s + "\n")
weather_path.each_child do |child|
  input.write((child.dirname.dirname.basename + child.dirname.basename + child.basename).to_s + "\n")
end
