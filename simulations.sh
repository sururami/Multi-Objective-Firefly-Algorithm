rm -rf 1ZDD1/
rm -rf prot.*

cd src/ && make &&  cd ..  &&  ./Mo_Fa_main CHARMm22 4 2 34 1000 8 50000 instances/1ZDD.seq 1ZDD1 1
