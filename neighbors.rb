offsets = 
"0 0 0
1 1 0
-1 -1 0
1 -1 0
-1 1 0
1 0 1
-1 0 -1
1 0 -1
-1 0 1
0 1 1
0 -1 -1
0 1 -1
0 -1 1".split("\n").map { |s| s.split.map(&:to_i) }

# find all next nearest neighbors of 0, 0, 0
nnn = offsets.dup.map! { |node| offsets.map { |o| [node, o].transpose.map { |a| a.inject(:+) } } }.flatten!(1).uniq! - offsets;

# find unique parallel vectors
nnn.each do |a|
  nnn.delete_if { |b| a == b.map(&:-@)}
end

nnn.sort!.reverse!

puts nnn.map { |a| a.flatten.join(" ") }
