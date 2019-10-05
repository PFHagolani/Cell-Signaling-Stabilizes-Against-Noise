#!/bin/bash


cd src/core

pwd

aleas=aleas.mod.f90

#making object files
echo 'making object files'

gfortran -S -w OpenGL_gl.f90 OpenGL_glu.f90 OpenGL_glut.f90
gfortran  -c $aleas
gfortran  -c general.mod.f90 
gfortran  -c geompack3.f90 gnuplotter.mod.f90
gfortran  -c shell.mod.f90 neighboring.mod.f90
gfortran  -c genetic.mod.f90
gfortran  -c energy.mod.f90 io.mod.f90 
gfortran -w -c biomechanic.mod.f90 death.mod.f90 pola.mod.f90 ecm.mod.f90 growth.mod.f90
gfortran  -c ic.mod.f90 mitosis.mod.f90
gfortran  -c single_node.mod.f90
gfortran  -c conservative6.mod.f90
gfortran  -c complexity6.mod.f90
gfortran  -c fitmo.mod.f90
gfortran  -c nexus.mod.f03
gfortran  -c inicial.mod.f90
gfortran -w -c editor.mod.f90 model.mod.f90
gfortran -w -c drawer.mod.f90 automaticon.mod.f90
gfortran  -c elli.f90
 

#linking
echo 'linking'

gfortran -w gnuplotter.mod.o aleas.mod.o general.mod.o neighboring.mod.f90 genetic.mod.o energy.mod.o shell.mod.o io.mod.o pola.mod.o mitosis.mod.o growth.mod.o death.mod.o single_node.mod.o ic.mod.o ecm.mod.o complexity6.mod.o conservative6.mod.o fitmo.mod.o nexus.mod.o biomechanic.mod.o model.mod.o inicial.mod.o editor.mod.o drawer.mod.f90 automaticon.mod.f90 geompack3.f90 elli.o -o EMaker -lGL -lGLU -lglut

#cleaning
echo 'cleaning'
rm *.s *.o *.mod

cd ../..
mv src/core/EMaker bin 



