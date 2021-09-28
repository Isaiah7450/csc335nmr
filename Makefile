# Created by: Isaiah Hoffman
# Created on: October 20, 2020
.RECIPEPREFIX := >
.PHONY := clean rebuild
Prog_Name := root_finder

My_Flags := -g -Wall -Wextra -pedantic-errors

$(Prog_Name).out: $(Prog_Name).o main.o
> gfortran $(My_Flags) $(Prog_Name).o main.o -o $(Prog_Name).out

$(Prog_Name).o: $(Prog_Name).f95
> gfortran $(My_Flags) -c $(Prog_Name).f95 -o $(Prog_Name).o

main.o: main.f95
> gfortran $(My_Flags) -c main.f95 -o main.o

clean:
> rm -f *.o
> rm -f *.mod
> rm -f *.out

rebuild: clean $(Prog_Name).out

