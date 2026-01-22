default:
	@echo 'No default for this Makefile, run make help for targets'
	
help:
	@echo "------------------------------------------------------------------------------------"
	@echo "Targets:"
	@echo "   clean : cleans basilisk"
	@echo "   download : downloads basilisk from basilisk website using darcs"
	@echo "   update : updates basilisk from basilisk website using darcs"
	@echo "   compile : compiles basilisk with the selected architecture"
	@echo "------------------------------------------------------------------------------------"

 
clean:
ifeq ($(shell test -d $(BASILISK) && echo -n yes),yes)
	$(RM) $(BASILISK_COMPILED)
	@echo "Compiled basilisk $(BASILISK_COMPILED) deleted"
else
	@echo "No compiled basilisk $(BASILISK_COMPILED) to delete"
endif


download:
ifeq ($(shell test -d $(BASILISK_HEAD) && echo -n yes),yes)
	@echo "Basilisk head directory $(BASILISK_HEAD) already exists, do \"make update\" instead"
else
	darcs get --lazy http://basilisk.fr/basilisk
endif

	
update:
ifeq ($(shell test -d $(BASILISK_HEAD) && echo -n yes),yes)
	cd $(BASILISK_HEAD); darcs pull; cd $(OCTREE_HOME)
else
	@echo "Basilisk head directory $(BASILISK_HEAD) does not exist yet, do \"make download\" first"
endif


compile: clean	
	$(MKDIR) $(BASILISK_COMPILED)
	$(CPDIR) $(BASILISK_HEAD)/src $(BASILISK_COMPILED)/
	cd $(BASILISK); $(CPDIR) config.gcc config; make -k clean; make; cd $(OCTREE_HOME)
	
		
##################################################################
# internal commands                                              #
##################################################################
RM    := rm -fr
MKDIR := mkdir -p
EXPAND := tar xzvf
CPDIR := cp -r

