#!/bin/sh
for T in 0.034 0.035 0.036 0.037 0.038
do
    name=kmc-default-T$T
    mkdir "$name"
    cp kmc.o $name/
    pushd $name
    echo "#!/bin/sh" > sub
    echo "#$ -S /bin/sh" >> sub
    echo "#$ -N pmf_0.5_0.0_1.0" >> sub
    echo "#$ -cwd" >> sub
    echo "#$ -j y" >> sub
    echo "./kmc.o -o $name -T $T" >> sub
    qsub sub
    popd
done
