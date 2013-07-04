# This is a script written 2013-07-03 to determine if the xml file is providing
# Reasonable or unreasonable outputs.

require 'fileutils'
require 'pathname'
require 'date'

if __FILE__ == $0
  filename = "output.xml"
  filename = ARGV.first if ARGV.first
  print "Input File: ", filename, "\n"
  f = File.new(filename, 'r')
  f.each_line.with_index do |line, line_number|
    if line.split[2]
      word = line.split[2]
      probability = word.split("=")
      if probability[0] == "probability"
        prob_number = probability[1].gsub(/[^0-9A-Za-z\.]/, '').to_f
        print line_number, " ", line.strip, "\n" if prob_number < 0.2
      end
    end
  end
end
