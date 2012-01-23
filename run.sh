#!/bin/sh
# ./kmc.o -T 0.037757 -T 0.0375 -t 10000 -tp 10000 -n 10000000 -i 10000 -dm 3000 -dpm 3000 -o test.xyz # 10k
# ./kmc.o  -T 0.0350 -t 10000 -tp 10000 -n 10000000 -i 10000 -dm 3000 -dpm 3000 -o test.xyz -e efile-michael # michael's scale
./kmc.o  -T 0.0350 -t 2000 -tp 200 -n 1000000 -i 1000 -dm 3000 -dpm 3000 -o test.xyz -e efile-michael # comp new
# ./kmc.o -T 0.0370 -t 10000 -tp 0 -n 10000000 -i 10000 -dm 3000 -dpm 0 -o test.xyz # oct
# ./kmc.o -T 0.0345 -t 10000 -tp 0 -n 10000000 -i 10000 -dm 3000 -dpm 0 -o test.xyz #truncated oct
# ./kmc.o -T 0.0375 -t 2000 -tp 2000 -n 10000000 -i 10000 -o test.xyz # low number
# ./kmc.o -T 0.0375 -t 2000 -tp 2000 -n 10000000 -i 10000 -o test.xyz # low number low pvp
# ./kmc.o -T 0.0385 -t 10000 -tp 0 -n 1000000 -i 1000 -dm 3000 -dpm 0 -o test.xyz # no pvp
