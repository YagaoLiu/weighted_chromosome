Compile with `make`

run pattern matching with 
```console
 (/usr/bin/time -v ./index -t ../chr22.in -z 2 -p ../z2 -o pm_time/z2time.txt) &> peak_ram/z2peak.txt
 (/usr/bin/time -v ./index -t ../chr22.in -z 4 -p ../z4 -o pm_time/z4time.txt) &> peak_ram/z4peak.txt
 (/usr/bin/time -v ./index -t ../chr22.in -z 8 -p ../z8 -o pm_time/z8time.txt) &> peak_ram/z8peak.txt
 (/usr/bin/time -v ./index -t ../chr22.in -z 16 -p ../z16 -o pm_time/z16time.txt) &> peak_ram/z16peak.txt
 (/usr/bin/time -v ./index -t ../chr22.in -z 32 -p ../z32 -o pm_time/z32time.txt) &> peak_ram/z32peak.txt
```
try Index size:
```console 
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 2 -p ../z2
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 4 -p ../z4
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 8 -p ../z8
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 16 -p ../z16
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 32 -p ../z32
```
