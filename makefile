mvmc:
	cd src/mVMC/;make -f makefile_src 
	cd tool;make -f makefile_tool

help:
	@echo ""
	@echo "Usage :"
	@echo "make <entry>"
	@echo ""
	@echo "<entry> is chosen from below"
	@echo "      mvmc : Build simulator mVMC in src/ and tool/"
	@echo " userguide : Generate userguid_jp.pdf & userguide_en.pdf in doc/"
	@echo "     clean : Remove all generated files excepting makefile and doc/"
	@echo " veryclean : Remove all generated files including makefile and doc/"
	@echo ""

userguide:
	cd doc/jp; make html latexpdfja
	cd doc/en; make html latexpdfja

new:
	cd src/mVMC/; rm -f *.o vmc.out vmcdry.out
	cd ../../;
	make
	make install 

luna:
	cp src/make.sys.luna src/make.sys
	make new

sekirei:
	cp src/make.sys.sekirei src/make.sys
	make new

clean:
	cd src/mVMC/; make -f makefile_src clean
	cd tool; make -f makefile_tool clean

veryclean:
	make clean
	cd doc/jp; make clean
	cd doc/en; make clean
	rm -f src/make.sys makefile

install:
	cp src/mVMC/vmc.out ~/bin/dvmc.out
	cp src/mVMC/vmcdry.out ~/bin/dvmcdry.out

