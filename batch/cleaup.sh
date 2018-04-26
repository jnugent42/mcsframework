#!/bin/bash

find batchscript_* -not -name batchscript_XX -exec rm {} \;
#rm job*
rm *.out
rm *.err
rm *.log
find sub_script_* -not -name sub_script_XX.sh -exec rm {} \;
find mcbatchscript_* -not -name mcbatchscript_YY -exec rm {} \;
find mcsub_script_* -not -name mcsub_script_YY.sh -exec rm {} \;
