if [[ -z "$KINFIT2PATH" ]]
then
    echo "Trying to initialize environment:"
    source ./setup.sh
fi
if [[ -z "$KINFIT2PATH" ]]
then
    echo "Trying to initialize environment failed!"
else
    echo "removing old files"
    rm -f libHHKinFit2Core.so
    rm -f libHHKinFit2Scenarios.so
    
    echo "creating shared library"
    g++ -Wall -fPIC -shared HHKinFit2Core/src/*.cpp `root-config --cflags --glibs` -I $KINFIT2PATH -o libHHKinFit2Core.so 
    g++ -Wall -fPIC -shared HHKinFit2Scenarios/src/*.cpp `root-config --cflags --glibs` -I $KINFIT2PATH -o libHHKinFit2Scenarios.so 

    echo "creating executable"
    #g++ -std=c++11 bin/compare2DKinFits.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2/HHKinFit2 -lHHKinFit2 -L $KINFIT1PATH -lHHKinFit1 -o compare2DKinFits
    g++ HHKinFit2Testing/bin/fitSingleEvent2DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2 -lHHKinFit2Core -lHHKinFit2Scenarios -o fitSingleEvent2DKinFit
    g++ HHKinFit2Testing/bin/fitSingleEvent1DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH -L $KINFIT2PATH/HHKinFit2 -lHHKinFit2Core -lHHKinFit2Scenarios -o fitSingleEvent1DKinFit

    # clang version -> paths not correct!!
    # clang -Wall -Wextra -fPIC -shared src/*.cpp `root-config --cflags --glibs` -I ./interface -D HHKINFIT2 -o libHHKinFit2.so 
    # clang examples/fitSingleEvent2DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH/interface -L $KINFIT2PATH -lHHKinFit2 -stdlib=libstdc++ -lstdc++ -D HHKINFIT2 -o fitSingleEvent2DKinFit
    # clang examples/fitSingleEvent1DKinFit.C `root-config --cflags --glibs` -I $KINFIT2PATH/interface -L $KINFIT2PATH -lHHKinFit2 -stdlib=libstdc++ -lstdc++ -D HHKINFIT2 -o fitSingleEvent1DKinFit
fi


