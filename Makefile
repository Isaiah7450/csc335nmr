# Created by: Isaiah Hoffman
# Created on: October 20, 2020
.RECIPEPREFIX := >
.PHONY := clean rebuild
Prog_Name := nmr

My_Flags := -g -Wall -Wextra -pedantic-errors

$(Prog_Name).out: type_library.o root_finder.o interpolation.o \
  calculus.o main.o
> gfortran $(My_Flags) type_library.o root_finder.o main.o -o $(Prog_Name).out

type_library.o: type_library.f95
> gfortran $(My_Flags) -c type_library.f95 -o type_library.o

root_finder.o: root_finder.f95 type_library.o
> gfortran $(My_Flags) -c root_finder.f95 -o root_finder.o

interpolation.o: interpolate.f95 type_library.o
> gfortran $(My_Flags) -c interpolate.f95 -o interpolation.o

calculus.o: calculus.f95 type_library.o
> gfortran $(My_Flags) -c calculus.f95 -o calculus.o

main.o: main.f95 root_finder.o interpolation.o type_library.o
> gfortran $(My_Flags) -c main.f95 -o main.o

clean:
> rm -f *.o
> rm -f *.mod
> rm -f *.out

rebuild: clean $(Prog_Name).out

