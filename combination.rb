# ruby file to generate all of the permutation arrays

def main()
  a = (0..8).to_a
  
  k = 4
  count = 0
  1.upto(k) do |i|
    enum = a.combination(i)
    
    enum.each do |array|
      current_array = array
      while current_array.length < k
        current_array << -1
      end
      count += 1
      
      print "{"
      current_array.each do |elt|
        print elt.to_s + ","
      end
      print "}, "
    end
  end
  print "\n\n" + count.to_s + "\n"
end

main()