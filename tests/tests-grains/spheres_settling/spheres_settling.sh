#!/bin/bash
./grains_seq Grains/Init/insert.xml 1> /dev/null 2> log
MAX_VEL=$(awk '{print $2}' Grains/Init/insert_VitesseMaxMean.dat | tail -1)
SUCCESS=$(python3 -c $"if ${MAX_VEL}<1.e-10:"$'\n'$"    print(\"0\")"$'\n'$"else:"$'\n'$"    print(\"1\")")
echo "$SUCCESS" | grep -q "0"
