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
  weather_dir = "RandomWeather"

  ARGV.each_with_index do |arg, index|
    seed                = ARGV[index+1].to_i if arg == "-seed"
    number_of_members   = ARGV[index+1].to_i if arg == "-num_members"
    num_weather_points  = ARGV[index+1].to_i if arg == "-num_weather_points"
    max_lat             = ARGV[index+1].to_f if arg == "-max_lat"
    min_lat             = ARGV[index+1].to_f if arg == "-min_lat"
    max_prob            = ARGV[index+1].to_f if arg == "-max_prob"
    min_prob            = ARGV[index+1].to_f if arg == "-min_prob"
    weather_dir         = ARGV[index+1] if arg == "-weather_dir"
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
      f.write("# Product type: CONVECTION\n")
      f.write("# Variable: DevProb\n")
      f.write("# Unit: scalar\n")
      f.write("# Input filename: null\n")
      f.write("# Member: 1/6\n")
      f.write("# Probability: ")
      if member_number < number_of_members
        f.write(file_probabilities[member_number].to_s) 
        probability_written += file_probabilities[member_number]
      else if member_number == number_of_members
        f.write((1 - probability_written).to_s)
      end
      f.write("\n# Reference time: 20090618T060000\n")
      f.write("# Valid time, min valid time range, max valid time range: 20090618T060000,20090618T060000,20090618T061500\n")
      f.write("# Grid size x,y,z: 322, 273, 5\n")
      f.write("# Grid spacing (meters,meters,ft): 5000.00000,5000.00000,5000.00000\n")
      f.write("# Dataset: x,y,z,lat,lon,alt,value\n")
      
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
end
