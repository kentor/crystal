#!/bin/sh
./kmc_p.o -T 0.0375 -t 10000 -tp 10000 -n 10000000 -i 10000 -dm 3000 -dpm 3000 -o test.xyz # 10k
# ./kmc_p.o -T 0.0375 -t 2000 -tp 2000 -n 10000000 -i 10000 -o test.xyz # low number
# ./kmc_p.o -T 0.0375 -t 2000 -tp 2000 -n 10000000 -i 10000 -o test.xyz # low number low pvp
# ./kmc_p.o -T 0.0375 -t 10000 -tp 0 -n 10000000 -i 10000 -dm 3000 -dpm 0 -o test.xyz # no pvp
