#!/bin/tcsh
foreach x(`seq 1 1 500`)   
    setenv iCheck `more log.txt | grep "Zee_53X_SIM_PAT_"$x"_" | wc -l`
    if( $iCheck == "2" ) then
        echo $x
    endif
end
unsetenv iCheck
