valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 2 -p ../z2
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 4 -p ../z4
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 8 -p ../z8
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 16 -p ../z16
valgrind --tool=massif --time-unit=ms --massif-out-file=massif.out.wst.z2 ./index -t ../chr22.in -z 32 -p ../z32
