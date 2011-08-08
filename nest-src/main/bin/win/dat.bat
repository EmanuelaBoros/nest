@echo off

set NEST_HOME=%CD%

IF ["%NEST_HOME%"]==[] echo "NEST_HOME is not defined. Please set NEST_HOME=nest_installation_folder"

IF [%NEST_HOME:~-1%]==[/] set NEST_HOME=%NEST_HOME:~0,-1%
IF [%NEST_HOME:~-1%]==[\] set NEST_HOME=%NEST_HOME:~0,-1%

"%NEST_HOME%\jre\bin\java.exe" ^
    -server -Xmx850M -XX:CompileThreshold=100 -Xverify:none ^
    -XX:+AggressiveOpts -XX:+UseFastAccessorMethods -Xconcurrentio ^
    -Dceres.context=nest ^
    -Dnest.debug=false ^
    -Dnest.consoleLog=true ^
    -Dnest.logLevel=INFO ^
    "-Dnest.home=%NEST_HOME%" ^
    -jar "%NEST_HOME%\bin\ceres-launcher.jar" -d %*

exit /B 0
