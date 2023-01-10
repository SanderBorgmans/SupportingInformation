


for i in *fold/; do
    cd $i
    for j in 1 2 4 8 16 32 64; do
        if [ ! -d "loading_${j}" ]; then mkdir loading_${j}; fi
        
        cd loading_${j}
        cp ../../H2O.chk .
        cp ../../*.py .
        cp ../../pars*.txt .
        cp ../framework.chk .

        python3 insert_molecules.py $j
        cd ../
    done
    cd ../
done

