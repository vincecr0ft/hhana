
# student name
HHSTUDENT ?= hhskim
# ntuple production directory
HHNTUP ?= ntuples/prod_v29/hhskim
# ntuple running directory
HHNTUP_RUNNING ?= /afs/cern.ch/user/v/vcroft/portal/hhskims/qskims
# maximum number of processors to request in PBS
PBS_PPN_MAX ?= 15

# current git branch
BRANCH := $(shell git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')

.PHONY: dump

default: clean

clean-root:
	rm -f $(HHNTUP)/$(HHSTUDENT).root

clean-h5:
	rm -f $(HHNTUP)/$(HHSTUDENT).h5

clean-ntup: clean-root clean-h5

clean-grl:
	rm -f $(HHNTUP)/observed_grl_11.xml
	rm -f $(HHNTUP)/observed_grl_12.xml
	rm -f ~/observed_grl_11.xml
	rm -f ~/observed_grl_12.xml
	rm -f $(HHNTUP)/merged_grl_11.xml
	rm -f $(HHNTUP)/merged_grl_12.xml

clean-pyc:                                                                      
	find mva statstools -name "*.pyc" -exec rm {} \;

clean: clean-pyc

check-files:
	./checkfile $(HHNTUP_RUNNING)/$(HHSTUDENT)*.root

check-ntup:
	./checkfile $(HHNTUP)/$(HHSTUDENT).root

browse:
	rootpy browse $(HHNTUP)/$(HHSTUDENT).root

roosh:
	roosh $(HHNTUP)/$(HHSTUDENT).root

$(HHNTUP_RUNNING)/$(HHSTUDENT).data12-JetTauEtmiss.root:
	@if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).data12-JetTauEtmiss_1.root ]; then \
		test -d $(HHNTUP_RUNNING)/data || mkdir $(HHNTUP_RUNNING)/data; \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).data12-JetTauEtmiss_*.root $(HHNTUP_RUNNING)/data; \
		hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).data12-JetTauEtmiss.root $(HHNTUP_RUNNING)/data/$(HHSTUDENT).data12-JetTauEtmiss_*.root; \
		test -d $(HHNTUP_RUNNING)/data_log || mkdir $(HHNTUP_RUNNING)/data_log; \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).data12_*.e[0-9]* $(HHNTUP_RUNNING)/data_log/; \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).data12_*.o[0-9]* $(HHNTUP_RUNNING)/data_log/; \
		mv $(HHNTUP_RUNNING)/supervisor-$(HHSTUDENT)-$(HHSTUDENT).data12-JetTauEtmiss_*.log $(HHNTUP_RUNNING)/data_log/; \
	fi

$(HHNTUP_RUNNING)/$(HHSTUDENT).data11-JetTauEtmiss.root:
	@if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).data11-JetTauEtmiss_1.root ]; then \
		test -d $(HHNTUP_RUNNING)/data || mkdir $(HHNTUP_RUNNING)/data; \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).data11-JetTauEtmiss_*.root $(HHNTUP_RUNNING)/data; \
		hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).data11-JetTauEtmiss.root $(HHNTUP_RUNNING)/data/$(HHSTUDENT).data11-JetTauEtmiss_*.root; \
		test -d $(HHNTUP_RUNNING)/data_log || mkdir $(HHNTUP_RUNNING)/data_log; \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).data11_*.e[0-9]* $(HHNTUP_RUNNING)/data_log/; \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).data11_*.o[0-9]* $(HHNTUP_RUNNING)/data_log/; \
		mv $(HHNTUP_RUNNING)/supervisor-$(HHSTUDENT)-$(HHSTUDENT).data11-JetTauEtmiss_*.log $(HHNTUP_RUNNING)/data_log/; \
	fi

init-data-12: $(HHNTUP_RUNNING)/$(HHSTUDENT).data12-JetTauEtmiss.root

init-data-11: $(HHNTUP_RUNNING)/$(HHSTUDENT).data11-JetTauEtmiss.root

init-data: init-data-11 init-data-12

$(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_TES_FAKE_TOTAL_UP.root:
	@test -d $(HHNTUP_RUNNING)/embed_tes || mkdir $(HHNTUP_RUNNING)/embed_tes
	@for TES_TERM in TES_TRUE_MODELING TES_TRUE_SINGLEPARTICLEINTERPOL TES_TRUE_INSITUINTERPOL TES_FAKE_TOTAL TES_TRUE_TOTAL; do \
		if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_UP_1.root ]; then \
			mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_*.root $(HHNTUP_RUNNING)/embed_tes; \
		fi; \
		if [ -f $(HHNTUP_RUNNING)/embed_tes/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_UP_1.root ]; then \
			hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_UP.root $(HHNTUP_RUNNING)/embed_tes/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_UP_*.root; \
			hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_DOWN.root $(HHNTUP_RUNNING)/embed_tes/$(HHSTUDENT).embed12-HH-IM_$${TES_TERM}_DOWN_*.root; \
		fi; \
	done

$(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM.root:
	@test -d $(HHNTUP_RUNNING)/embed || mkdir $(HHNTUP_RUNNING)/embed
	@if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_1.root ]; then \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_[0-9]*.root $(HHNTUP_RUNNING)/embed; \
	fi
	@if [ -f $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed12-HH-IM_1.root ]; then \
		hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM.root $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed12-HH-IM_*.root; \
	fi
	@if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-UP_1.root ]; then \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-UP_*.root $(HHNTUP_RUNNING)/embed; \
	fi
	@if [ -f $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed12-HH-UP_1.root ]; then \
		hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-UP.root $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed12-HH-UP_*.root; \
	fi
	@if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-DN_1.root ]; then \
		mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-DN_*.root $(HHNTUP_RUNNING)/embed; \
	fi
	@if [ -f $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed12-HH-DN_1.root ]; then \
		hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-DN.root $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed12-HH-DN_*.root; \
	fi

$(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim_TES_TRUE_FINAL_UP.root:
	@test -d $(HHNTUP_RUNNING)/embed_tes || mkdir $(HHNTUP_RUNNING)/embed_tes
	@for TES_TERM in TES_TRUE_FINAL TES_FAKE_FINAL; do \
		if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_UP_1.root ]; then \
			mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_*.root $(HHNTUP_RUNNING)/embed_tes; \
		fi; \
		if [ -f $(HHNTUP_RUNNING)/embed_tes/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_UP_1.root ]; then \
			hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_UP.root $(HHNTUP_RUNNING)/embed_tes/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_UP_*.root; \
			hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_DOWN.root $(HHNTUP_RUNNING)/embed_tes/$(HHSTUDENT).embed11-hh-isol-mfsim_$${TES_TERM}_DOWN_*.root; \
		fi; \
	done

$(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim.root:
	@test -d $(HHNTUP_RUNNING)/embed || mkdir $(HHNTUP_RUNNING)/embed
	@for syst in im up dn; do \
		if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfs$${syst}_1.root ]; then \
			mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfs$${syst}_[0-9]*.root $(HHNTUP_RUNNING)/embed; \
		fi; \
		if [ -f $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed11-hh-isol-mfs$${syst}_1.root ]; then \
			hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfs$${syst}.root $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed11-hh-isol-mfs$${syst}_[0-9]*.root; \
		fi; \
	done
	@for syst in no tight; do \
		if [ -f $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-$${syst}isol-mfsim_1.root ]; then \
			mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-$${syst}isol-mfsim_[0-9]*.root $(HHNTUP_RUNNING)/embed; \
		fi; \
		if [ -f $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed11-hh-$${syst}isol-mfsim_1.root ]; then \
			hadd $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-$${syst}isol-mfsim.root $(HHNTUP_RUNNING)/embed/$(HHSTUDENT).embed11-hh-$${syst}isol-mfsim_[0-9]*.root; \
		fi; \
	done

.PHONY: embed-log
embed-log:
	@test -d $(HHNTUP_RUNNING)/embed_log || mkdir $(HHNTUP_RUNNING)/embed_log
	@-mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed1[1-2]-*.e[0-9]* $(HHNTUP_RUNNING)/embed_log/
	@-mv $(HHNTUP_RUNNING)/$(HHSTUDENT).embed1[1-2]-*.o[0-9]* $(HHNTUP_RUNNING)/embed_log/
	@-mv $(HHNTUP_RUNNING)/supervisor-$(HHSTUDENT)-$(HHSTUDENT).embed1[1-2]-*.log $(HHNTUP_RUNNING)/embed_log/

init-embed-11-nominal: $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim.root

init-embed-11-sys: $(HHNTUP_RUNNING)/$(HHSTUDENT).embed11-hh-isol-mfsim_TES_TRUE_FINAL_UP.root

init-embed-11: init-embed-11-nominal init-embed-11-sys

init-embed-12-sys: $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM_TES_FAKE_TOTAL_UP.root

init-embed-12-nominal: $(HHNTUP_RUNNING)/$(HHSTUDENT).embed12-HH-IM.root

init-embed-12: init-embed-12-nominal init-embed-12-sys embed-log

init-embed: init-embed-11 init-embed-12 embed-log

init-mc:
	@test -d $(HHNTUP_RUNNING)/mc_log || mkdir $(HHNTUP_RUNNING)/mc_log
	@-mv $(HHNTUP_RUNNING)/$(HHSTUDENT).*mc1[1-2]*.e[1-9]* $(HHNTUP_RUNNING)/mc_log/
	@-mv $(HHNTUP_RUNNING)/$(HHSTUDENT).*mc1[1-2]*.o[1-9]* $(HHNTUP_RUNNING)/mc_log/
	@-mv $(HHNTUP_RUNNING)/supervisor-$(HHSTUDENT)-$(HHSTUDENT).*mc1[1-2]*.log $(HHNTUP_RUNNING)/mc_log/

init-ntup: init-data init-embed init-mc

$(HHNTUP)/$(HHSTUDENT).root:
	@./merge-ntup -s $(HHSTUDENT) -o $(HHNTUP)/$(HHSTUDENT).root $(HHNTUP)/$(HHSTUDENT).*.root

$(HHNTUP)/$(HHSTUDENT).h5: $(HHNTUP)/$(HHSTUDENT).root
	@root2hdf5 --complib lzo --complevel 0 --quiet $^

ntup: $(HHNTUP)/$(HHSTUDENT).h5

.PHONY: ntup-update
ntup-update:
	@./merge-ntup -s $(HHSTUDENT) -o $(HHNTUP)/$(HHSTUDENT).root $(HHNTUP_RUNNING)/$(HHSTUDENT).*.root
	@root2hdf5 --update --complib lzo --complevel 0 --quiet $(HHNTUP)/$(HHSTUDENT).root

.PHONY: higgs-pt
higgs-pt:
	./higgs-pt $(HHNTUP_RUNNING)/hhskim*tautau*.root

.PHONY: $(HHNTUP)/merged_grl_11.xml
$(HHNTUP)/merged_grl_11.xml:
	ls $(HHNTUP)/data/hhskim.data11-*.root | sed 's/$$/:\/lumi/g' | xargs grl or > $@

.PHONY: $(HHNTUP)/merged_grl_12.xml
$(HHNTUP)/merged_grl_12.xml:
	ls $(HHNTUP)/data/hhskim.data12-*.root | sed 's/$$/:\/lumi/g' | xargs grl or > $@

.PHONY: $(HHNTUP)/observed_grl_11.xml
$(HHNTUP)/observed_grl_11.xml: $(HHNTUP)/merged_grl_11.xml ../higgstautau/grl/2011/current.xml 
	grl and $^ > $@

.PHONY: $(HHNTUP)/observed_grl_12.xml
$(HHNTUP)/observed_grl_12.xml: $(HHNTUP)/merged_grl_12.xml ../higgstautau/grl/2012/current.xml 
	grl and $^ > $@

.PHONY: ~/observed_grl_11.xml
~/observed_grl_11.xml: $(HHNTUP)/observed_grl_11.xml
	cp $^ $@

.PHONY: ~/observed_grl_12.xml
~/observed_grl_12.xml: $(HHNTUP)/observed_grl_12.xml
	cp $^ $@

.PHONY: grl-11
grl-11: ~/observed_grl_11.xml

.PHONY: grl-12
grl-12: ~/observed_grl_12.xml

grl: grl-11 grl-12

bundle:
	@rm -f ~/higgstautau-mva-plots.tar.gz
	@find plots/variables/$(BRANCH) -name '*.eps' -print0 | tar -vpcz --null -T - -f ~/higgstautau-mva-plots.tar.gz
	@echo bundle at ~/higgstautau-mva-plots.tar.gz

png-bundle:
	@rm -f ~/higgstautau-mva-plots.tar.gz
	@find plots/variables/$(BRANCH) -name '*.png' -print0 | tar -vpcz --null -T - -f ~/higgstautau-mva-plots.tar.gz
	@echo bundle at ~/higgstautau-mva-plots.tar.gz

montage:
	@montage -tile 4x5 -geometry 400x400+3+3 plots/*.png montage.pdf

thesis:
	@rm -f ~/thesis.tar.gz
	@find plots/contours/$(BRANCH) \
	      plots/categories/$(BRANCH) \
	      plots/fake_shapes/$(BRANCH) \
	      plots/shapes/$(BRANCH)/higgs_vs_ztautau \
	      plots/shapes/$(BRANCH)/qcd_vs_ztautau \
	      plots/normalization/$(BRANCH) \
	      plots/bdt/$(BRANCH) \
	      plots/variables/$(BRANCH) \
	      -maxdepth 1 -name '*.eps' -print0 | tar -vpcz --null -T - -f ~/thesis.tar.gz
	@echo created ~/thesis.tar.gz

workspace-plots:
	@rm -f ~/workspace.tar.gz
	@find ./workspaces -name "*.png" -print0 | tar -vpcz --null -T - -f ~/workspace.tar.gz
	@echo created ~/workspace.tar.gz

test:
	nosetests -s -v mva

.PHONY: norms
norms:
	@for year in 2011 2012; do \
		for model in OS_NONISOL nOS nOS_ISOL nOS_NONISOL SS SS_ISOL SS_NONISOL NONISOL; do \
			nohup ./norm --fakes-region $${model} --year $${year} & \
			nohup ./norm --no-embedding --fakes-region $${model} --year $${year} & \
		done; \
		nohup ./norm --fakes-region SS_NTRK --target-region OS_NTRK --year $${year} & \
		nohup ./norm --no-embedding --fakes-region SS_NTRK --target-region OS_NTRK --year $${year} & \
		nohup ./norm --fakes-region SS --target-region OS --year $${year} & \
		nohup ./norm --no-embedding --fakes-region SS --target-region OS --year $${year} & \
	done

.PHONY: model-plots
model-plots:
	@for year in 2011 2012; do \
		for model in nOS nOS_ISOL nOS_NONISOL SS SS_ISOL SS_NONISOL OS_NONISOL NONISOL; do \
			PBS_LOG=log run-cluster ./plot-features --fakes-region $${model} --year $${year} --categories presel; \
			for category in vbf boosted; do \
				PBS_LOG=log run-cluster ./plot-bdt --year $${year} --category-names $${category} --fakes-region $${model}; \
				PBS_LOG=log run-cluster ./plot-features --year $${year} --category-names $${category} --fakes-region $${model}; \
			done; \
		done; \
	done

.PHONY: plots
plots:
	@for year in 2011 2012; do \
		for category in vbf boosted rest; do \
			PBS_LOG=log PBS_MEM=12gb run-cluster ./plot-features --unblind --systematics --show-ratio --year $${year} --category-names $${category} --output-formats eps png; \
		done; \
		PBS_LOG=log PBS_MEM=12gb run-cluster ./plot-features --unblind --systematics --show-ratio --year $${year} --categories presel --output-formats eps png; \
	done

.PHONY: mva-plots
mva-plots:
	@for year in 2011 2012; do \
		for category in vbf boosted; do \
			PBS_LOG=log PBS_MEM=12gb run-cluster ./plot-bdt --unblind --systematics --year $${year} --category-names $${category} --output-formats eps png; \
		done; \
	done

.PHONY: train-vbf
train-vbf:
	@PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./train vbf --masses 125 --procs $(PBS_PPN_MAX) --max-fraction 0.15 --min-fraction-steps 100

.PHONY: train-boosted
train-boosted:
	@PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./train boosted --masses 125 --procs $(PBS_PPN_MAX) --learning-rate 0.01 --max-fraction 0.04 --min-fraction-steps 100

.PHONY: train
train: train-vbf train-boosted

.PHONY: train-vbf-each-mass
train-vbf-each-mass:
	@for mass in $$(seq 100 5 150); do \
		PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./train vbf --masses $${mass} --procs $(PBS_PPN_MAX); \
	done

.PHONY: train-boosted-each-mass
train-boosted-each-mass:
	@for mass in $$(seq 100 5 150); do \
		PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./train boosted --masses $${mass} --procs $(PBS_PPN_MAX); \
	done

.PHONY: train-each-mass
train-each-mass: train-vbf-each-mass train-boosted-each-mass

.PHONY: binning
binning:
	@for year in 2011 2012; do \
		for category in boosted vbf; do \
			PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./optimize-binning --year $${year} --categories $${category} --procs $(PBS_PPN_MAX); \
		done; \
	done

.PHONY: binning-each-mass
binning-each-mass:
	@for year in 2011 2012; do \
		for mass in $$(seq 100 5 150); do \
			PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./optimize-binning --year $${year} --mass $${mass} --procs $(PBS_PPN_MAX); \
		done; \
	done

.PHONY: cuts-workspaces
cuts-workspaces:
	@for mass in $$(seq 100 5 150); do \
		PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace cuts --systematics --unblind --years 2012 --categories cuts --masses $${mass}; \
		# PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace cuts --systematics --unblind --years 2011 --categories cuts_2011 --masses $${mass}; \
	done

.PHONY: cuts-sideband-workspace
cuts-sideband-workspace:
	@PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace cuts --sideband --systematics --unblind --years 2012 --categories cuts --masses 125; \

.PHONY: mva-workspaces
mva-workspaces:
	@for year in 2012 2011; do \
		for mass in $$(seq 100 5 150); do \
			PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace mva --systematics --unblind --years $${year} --masses $${mass} --clf-mass 125; \
		done; \
	done

.PHONY: mva-workspaces-binning-test
mva-workspaces-binning-test:
	PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace mva --systematics --unblind --years 2012 --masses 125 --clf-mass 125 --output-suffix default; \
	PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace mva --systematics --unblind --years 2012 --masses 125 --clf-mass 125 --clf-bins 10 --output-suffix 10bins; \
	PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace mva --systematics --unblind --years 2012 --masses 125 --clf-mass 125 --clf-bins 20 --output-suffix 20bins; \



.PHONY: mva-sideband-workspace
mva-sideband-workspace:
	@for year in 2012 2011; do \
		PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace mva --sideband --systematics --unblind --years $${year} --masses 125 --clf-mass 125; \
	done

.PHONY: mva-workspaces-multibdt
mva-workspaces-multibdt:
	@for year in 2012 2011; do \
		for mass in $$(seq 100 5 150); do \
			PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace mva --output-suffix multibdt --systematics --unblind --years $${year} --masses $${mass}; \
		done; \
	done

.PHONY: weighted-mass-workspaces
weighted-mass-workspaces:
	@for year in 2012 2011; do \
		for mass in 110 125 150; do \
			PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace weighted-mass --systematics --years $${year} --masses $${mass}; \
		done; \
	done

.PHONY: weighted-mass-workspaces-cba
weighted-mass-workspaces-cba:
	@for mass in 110 125 150; do \
		PBS_LOG=log PBS_MEM=18gb run-cluster ./workspace weighted-mass-cba --categories cuts --systematics --years 2012 --masses $${mass}; \
	done

.PHONY: workspaces
workspaces: mva-workspaces cuts-workspaces

.PHONY: combine-mva
combine-mva:
	@cd workspaces/hh_nos_nonisol_ebz_mva; \
	for mass in $$(seq 100 5 150); do \
		combine hh_11_vbf_$${mass} hh_12_vbf_$${mass} --name hh_vbf_$${mass}; \
		combine hh_11_boosted_$${mass} hh_12_boosted_$${mass} --name hh_boosted_$${mass}; \
		combine hh_11_combination_$${mass} hh_12_combination_$${mass} --name hh_combination_$${mass}; \
	done;

.PHONY: combine-cuts
combine-cuts:
	@cd workspaces/hh_nos_nonisol_ebz_cuts; \
	for mass in $$(seq 100 5 150); do \
		combine hh_12_cuts_boosted_loose_$${mass} hh_12_cuts_boosted_tight_$${mass} --name hh_12_cuts_boosted_$${mass}; \
		combine hh_12_cuts_vbf_lowdr_$${mass} hh_12_cuts_vbf_highdr_loose_$${mass} hh_12_cuts_vbf_highdr_tight_$${mass} --name hh_12_cuts_vbf_$${mass}; \
	done;

.PHONY: pruning
pruning:
	@for ana in mva; do \
		cd workspaces; \
		mkdir pruning_$${ana}; \
		cd pruning_$${ana}; \
		for mass in $$(seq 100 5 150); do \
			cp -r ../hh_nos_nonisol_ebz_$${ana}/hh_combination_$${mass} .; \
			cp ../hh_nos_nonisol_ebz_$${ana}/hh_combination_$${mass}.root .; \
		done; \
		fix-workspace --quiet --suffix split_norm_shape hh_combination_[0-9][0-9][0-9]; \
		fix-workspace --quiet --suffix drop_others_shapes --drop-others-shapes hh_combination_[0-9][0-9][0-9]; \
		fix-workspace --quiet --suffix drop_others_shapes_prune_norms --drop-others-shapes --prune-norms hh_combination_[0-9][0-9][0-9]; \
		for thresh in 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.96 0.97 0.98 0.99 1; do \
			fix-workspace --quiet --drop-others-shapes --prune-norms --prune-shapes --chi2-thresh $${thresh} --suffix chi2_$${thresh} hh_combination_[0-9][0-9][0-9]; \
			fix-workspace --quiet --drop-others-shapes --prune-norms --symmetrize --prune-shapes --chi2-thresh $${thresh} --suffix chi2_$${thresh}_sym hh_combination_[0-9][0-9][0-9]; \
			fix-workspace --quiet --drop-others-shapes --prune-norms --symmetrize-partial --prune-shapes --chi2-thresh $${thresh} --suffix chi2_$${thresh}_part_sym hh_combination_[0-9][0-9][0-9]; \
		done; \
		cd ..; \
	done

.PHONY: fix-mva
fix-mva:
	# IMPORTANT: update pruning chi2 threshold from plots made from pruning routine above
	@PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./fix-workspace --quiet --symmetrize --prune-shapes --chi2-thresh 0.9 --drop-others-shapes --prune-norms workspaces/hh_nos_nonisol_ebz_mva

.PHONY: fix-cuts
fix-cuts:
	@PBS_LOG=log PBS_PPN=$(PBS_PPN_MAX) run-cluster ./fix-workspace --quiet --symmetrize --prune-shapes --chi2-thresh 0.9 --drop-others-shapes --prune-norms --prune-samples workspaces/hh_nos_nonisol_ebz_cuts


grid-scores:
	./plot-grid-scores --clf bdts/clf_boosted_125_NONISOL_ebz_12_0.pickle --name Boosted --grid bdts/clf_boosted_125_NONISOL_ebz_12_0_grid_scores.pickle --format eps
	./plot-grid-scores --clf bdts/clf_vbf_125_NONISOL_ebz_12_0.pickle --name VBF --grid bdts/clf_vbf_125_NONISOL_ebz_12_0_grid_scores.pickle --format eps
