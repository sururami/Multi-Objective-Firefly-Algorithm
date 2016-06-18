rm -rf 1ZDD1/
rm -rf prot.*

cd src/ && make &&  cd ..  &&  ./i-paes2 CHARMm27 4 2 34 1000 1 50000 instances/1ZDD.seq 1ZDD1 1
