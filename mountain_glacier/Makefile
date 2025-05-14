EXEC := exec_grad

F90 	:= gfortran
CC 		:= gcc
ADFirstAidKit:= ../ADFirstAidKit

$(EXEC): adStack.o adBinomial.o global_variables_b.o global_variables_d.o forward_b.o forward_d.o driver.o
		$(F90) -o $@ $^

%.o: %.f90 
		$(F90) -c $<

adBinomial.o:
		$(CC) -c $(ADFirstAidKit)/adBinomial.c

adStack.o : 
		$(CC) -c $(ADFirstAidKit)/adStack.c

forward_b.f90 global_variables_b.f90: forward.f90 global_variables.f90
		tapenade -reverse -adjmodulename %_adj -head "forward_problem(V)/(M)" forward.f90 global_variables.f90
forward_d.f90 global_variables_d.f90: forward.f90 global_variables.f90
		tapenade -tangent -tgtmodulename %_tgt -head "forward_problem(V)/(M)" forward.f90 global_variables.f90

.PHONY: clean
clean:
	$(RM) $(EXEC) *.o *.msg *.mod *_db.f90 *_b.f90 *_d.f90 *~ tapStack* *.txt

