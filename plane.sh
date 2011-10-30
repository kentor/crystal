#!/bin/sh
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7
do
    name=p3b$i
    mkdir "$name"
    cp simppvp.o simppvp.py $name/
    pushd $name
    python simppvp.py $i $name
    popd
done
