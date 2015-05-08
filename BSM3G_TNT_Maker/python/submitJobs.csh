counter=1
for l in `cat dataFiles.txt`; do
  cp miniAOD_runData_COPY.py miniAOD_runData.py 
  sed -i s:DATAFILE:$l:g miniAOD_runData.py
  sed -i s:OUTPUTFILE:"OutputTree_"$counter:g miniAOD_runData.py
  cmsRun miniAOD_runData.py
  counter=`expr $counter + 1`
  ls -lhrt
done
