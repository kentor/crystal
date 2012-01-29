Offsets = 
"1 1 0
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
0 -1 1".split("\n").map! { |l| l.split.map(&:to_i) }

p Offsets

# find all next nearest neighbors of 0, 0, 0
nnn = Offsets.dup.map! do |neigh|
  Offsets.map do |o|
    [neigh, o].transpose.map { |a| a.inject(:+) }
  end
end.flatten!(1).uniq! - Offsets;

# find unique parallel vectors
nnn.each do |v|
  nnn.delete_if { |u| u == v.map(&:-@) }
end

nnn.sort!.reverse!

p nnn

puts nnn.map! { |a| a.flatten.join(" ") }
