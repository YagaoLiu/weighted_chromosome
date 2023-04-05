# weighted_chromosome

unzip the chr22.in.gz file
```console
gzip -dk file.gz
```

compile using `make`

Run the following scripts to get patterns with cartain z:
```console
./index -t chr22.in -z 32 -o ./patterns/z32
./index -t chr22.in -z 64 -o ./patterns/z64
./index -t chr22.in -z 128 -o ./patterns/z128
./index -t chr22.in -z 256 -o ./patterns/z256
./index -t chr22.in -z 512 -o ./patterns/z512
./index -t chr22.in -z 1024 -o ./patterns/z1024
./index -t chr22.in -z 2048 -o ./patterns/z2048
```

The length of chr22, after removing the multiply 'N' at the beginning and ending, is 35194566. (35.2*10^6)



