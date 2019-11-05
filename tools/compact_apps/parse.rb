#!/usr/bin/ruby

a=`ls dump/solution*`.split("\n")
rank=Array.new

a.each do |x|
	rank.push(x.split(/solution\./)[1].to_i)
	
end

#puts rank

concurrency=rank.max+1

sol=Array.new

for i in 0..rank.max
	file = File.new("dump/solution.#{i}","r")
	while (b=file.gets)
		sol.push(b.to_f)
	end
	file.close
end

#puts sol

file = File.new("output/data.#{concurrency}","w")
for i in 0..sol.length-1
	file.printf("%-5d%-10f\n",i.to_f,sol[i])
end
file.close



