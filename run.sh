#!/bin/sh
./kmc.o -T 0.0385 -t 10000 -tp 5000 -n 10000000 -i 10000 -dm 3000 -dpm 3000 -o test.xyz # 10k
# ./kmc.o -T 0.0370 -t 10000 -tp 0 -n 10000000 -i 10000 -dm 3000 -dpm 0 -o test.xyz # oct
# ./kmc.o -T 0.0375 -t 2000 -tp 2000 -n 10000000 -i 10000 -o test.xyz # low number
# ./kmc.o -T 0.0375 -t 2000 -tp 2000 -n 10000000 -i 10000 -o test.xyz # low number low pvp
# ./kmc.o -T 0.0385 -t 10000 -tp 0 -n 10000000 -i 10000 -dm 3000 -dpm 0 -o test.xyz # no pvp
