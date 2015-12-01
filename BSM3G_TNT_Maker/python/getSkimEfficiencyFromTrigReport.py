#!/usr/bin/env python

import sys, re

from math import sqrt
def main(argv) :
    """
    Very simple script to get skimming efficiency
    from Trigger Report outputs
    
    Usage:
    ./getSkimEfficiencyFromLogs <mode> <FILELIST>
    mode: mu-tau (0/default) or e-tau (1) skimming scheme
    e.g.: ./getSkimEfficiencyFromLogs.py 1 LOGDIR/*.log'
    
    """

    # Get file list from argument vector
    alen = len(argv)
    if alen < 2:
        sys.exit()
    FLIST   = argv[1:]


    # mu-tau or e-tau filters
    ETAUSKIM = False
    EMUSKIM = False
    MUMUSKIM = False
    EESKIM = False
    TTSKIM = False
    TSKIM = False
    SUSYSKIM = False
    BSM3GTNTSKIM = False
    
    if int(argv[0]) == 1:
        ETAUSKIM = True
    elif int(argv[0]) == 2:
    	EMUSKIM = True
    elif int(argv[0]) == 3:
    	MUMUSKIM = True
    elif int(argv[0]) == 4:
    	EESKIM = True
    elif int(argv[0]) == 5:
    	TTSKIM = True
    elif int(argv[0]) == 6:
    	TSKIM = True
    elif int(argv[0]) == 7:
    	SUSYSKIM = True
    elif int(argv[0]) == 8:
    	BSM3GTNTSKIM = True

    ### Define list of filter modules here!!! ###
    PathName = 'p'
    Filters =  ['tauGenJets',
                'genTauDecaysToMuonCands','genTauDecaysToHadronsCands',
                'selectedGenTauDecaysToMuonEta','selectedGenTauDecaysToMuonPt',
                'selectedGenTauDecaysToHadronsEta','selectedGenTauDecaysToHadronsPt',
                'selectedPFTaus','selectedMuons','muTauPairs','selectedMuTauPairs']
    if ETAUSKIM:
        PathName = 'p'
        Filters =  ['tauGenJets',
                    'genTauDecaysToElectronCands','genTauDecaysToHadronsCands',
                    'selectedGenTauDecaysToElectronEta','selectedGenTauDecaysToElectronPt',
                    'selectedGenTauDecaysToHadronsEta','selectedGenTauDecaysToHadronsPt',
                    'selectedPFTaus','selectedElectrons','elecTauPairs','selectedElecTauPairs']
    elif EMUSKIM:
    	PathName = 'p'
	Filters = ['tauGenJets',
		   'genTauDecaysToMuonCands', 'genTauDecaysToElectronCands',
		   'selectedGenTauDecaysToMuonEta', 'selectedGenTauDecaysToMuonPt',
		   'selectedGenTauDecaysToElectronEta', 'selectedGenTauDecaysToElectronPt',
		   'selectedElectrons', 'selectedMuons', 'eMuPairs', 'selectedEMuPairs']
    elif MUMUSKIM:
    	PathName = 'p'
	Filters = ['tauGenJets',
		   'genTauDecaysToMuonCands', 'selectedGenTauDecaysToMuonEta',
		   'selectedGenTauDecaysToMuonPt', 'selectedGenTauMuonPairs',
		    'selectedMuons', 'muMuPairs', 'selectedmuMuPairs']
    elif EESKIM:
    	PathName = 'p'
	Filters = ['tauGenJets',
		   'genTauDecaysToElectronCands', 'selectedGenTauDecaysToElectronEta',
		   'selectedGenTauDecaysToElectronPt', 'selectedGenTauElectronPairs',
		   'selectedElectrons', 'elecElecPairs', 
		   'selectedelecElecPairs']
    elif TTSKIM:
    	PathName = 'p'
	Filters = ['selectedPFTaus', 'TauTauPairs', 'selectedTauTauPairs']
    elif TSKIM:
    	PathName = 'p'
	Filters = ['selectedLooseHPSPatTau', 'selectedTausInclusive']
    elif SUSYSKIM:
    	PathName = 'p'
	Filters = ['patConversions']
    elif BSM3GTNTSKIM:
    	PathName = 'p'
	Filters = ['TNT']
        

    ### Init counters ###
    NEvents = 0
    NPassedAll = 0


### Loop to get  total skim efficiency

    # Reg. Expr. patterns to find right lines
    fullEntryPattern = re.compile('TrigReport ---------- Event  Summary ------------[\t\n\r\f\v\W\w]+TrigReport ---------- Path   Summary')
    ppat = 'TrigReport Events total[\W\w]+'
    p = []
    for i in range(0,len(Filters)):
        p.append(re.compile(ppat))

    # Loop over all files
    for filen in FLIST:
        input = file(filen)
        fullTxt = input.read() # read full txt
        m0 = fullEntryPattern.search(fullTxt) # extract region
        if (m0):
            lines = re.split('\n',m0.group())
            for line in lines:
                for i in range(0,len(p)):
                    m = p[i].search(line)
                    if (m):
                        tabs = re.split('[\s]+',m.group())
                        if (i==0):
                            NEvents = NEvents + int(tabs[4])
                            NPassedAll = NPassedAll + int(tabs[7])
        input.close()

    exit



### Second loop on muCaloTauFilter

    ### Init counters ###
    NPassed = [0] * len(Filters)

    # Reg. Expr. patterns to find right lines
    fullEntryPattern = re.compile('TrigReport ---------- Modules in Path: '+PathName+' ----[\t\n\r\f\v\W\w]+TrigReport ---------- Module Summary')
    ppat = 'TrigReport [\W\w]+ '
    p = []
    count = []
    for i in range(0,len(Filters)):
        p.append(re.compile(ppat+Filters[i]))
        count.append(0)
        
    # Loop over all files
    for filen in FLIST:
        for i in range(0,len(Filters)):
            count[i] = 0
        input = file(filen)
        fullTxt = input.read() # read full txt
        m0 = fullEntryPattern.search(fullTxt) # extract region
        if (m0):
            lines = re.split('\n',m0.group())
            count[i] = 0
            for line in lines:
                for i in range(0,len(p)):
                    m = p[i].search(line)
                    if (m and count[i]==0):
                        count[i] = count[i] + 1
                        tabs = re.split('[\s]+',m.group())
                        NPassed[i] = NPassed[i] + int(tabs[4])            
        input.close()

    print '\n******************************'
    print PathName+' Efficiencies ***'
    print '******************************'
    print 'Events processed:',NEvents
    print 'Filter                              Passed   Efficiency   Cumul. Efficiency'
    print '-----------------------------------------------------------------------'
    eff = -1.
    # Treat first entry
    if (float(NEvents!=0)):
        eff = float(NPassed[0])/float(NEvents)
    print '%35s: %8d     %5.3f       %5.3f'%(Filters[0],NPassed[0],float(NPassed[0])/float(NEvents),float(NPassed[0])/float(NEvents))
    for i in range(1,len(Filters)):
        if (float(NPassed[i-1]==0)):
            eff = -1.
        else:
            eff = float(NPassed[i])/float(NPassed[i-1])                
        print '%35s: %8d     %5.3f       %5.3f'%(Filters[i],NPassed[i],eff,float(NPassed[i])/float(NEvents))
    print '-----------------------------------------------------------------------'
    last = len(Filters)-1
    print '%35s: %8d                 %5.3f'%('Total',NPassed[last],float(NPassed[last])/float(NEvents))
    
    newEff = float(NPassed[last])/float(NEvents)
    newEffError = sqrt(newEff*(1.0 - newEff)/float(NEvents))

    print '\n*** Eff / EffError: %8.6f %8.6f'%(newEff,newEffError)


    if ETAUSKIM:
        print '\n******************************'
        print ' Filter module descriptions ***'
        print '******************************'
        print '                       tauGenJets:  produce generator collection of taus '
        print '      genTauDecaysToElectronCands:  >= 1 tau which decays to electron '
        print '       genTauDecaysToHadronsCands:  >= 1 tau which decays hadronically '
        print 'selectedGenTauDecaysToElectronEta:  Electron from tau decay which has abs(eta)<2.5 '
        print ' selectedGenTauDecaysToElectronPt:  ... which has Pt>5. '
        print ' selectedGenTauDecaysToHadronsEta:  Jet from tau decay which has abs(eta)<2.5 '
        print '  selectedGenTauDecaysToHadronsPt:  ... which has Pt>10. '
        #print '                  selectedPFTaus:  >= 1 PFTau which pass fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut'
        print '                 selectedCaloTaus:  >= 1 CaloTau which pass '
        print '                                    caloRecoTauDiscriminationByLeadingTrackPtCut (lead trk pt>5.)'
        print '                selectedElectrons:  >= 1 pixelmatched gsf electron with'
        print '                                    pt > 10 & abs(eta) < 2.5 & eSuperClusterOverP>0.8'
        print '                                    & eSuperClusterOverP<1.25'
        #print '                  elecPFTauPairs:  produce elec-PFTau pairs '
        #print '      selectedElectronPFTauPairs:  require elec-PFTau pair separated by R>0.7 '
        print '                 elecCaloTauPairs:  produce elec-CaloTau pairs '
        print '         selectedElecCaloTauPairs:  require elec-CaloTau pair separated by R>0.7 '
        print '-----------------------------------------------------------------------'        
    else:
        print '\n******************************'
        print ' Filter module descriptions ***'
        print '******************************'
        print '                      tauGenJets:  produce generator collection of taus '
        print '         genTauDecaysToMuonCands:  >= 1 tau which decays to muon '
        print '      genTauDecaysToHadronsCands:  >= 1 tau which decays hadronically '
        print '   selectedGenTauDecaysToMuonEta:  Muon from tau decay which has abs(eta)<2.5 '
        print '    selectedGenTauDecaysToMuonPt:  ... which has Pt>5. '
        print 'selectedGenTauDecaysToHadronsEta:  Jet from tau decay which has abs(eta)<2.5 '
        print ' selectedGenTauDecaysToHadronsPt:  ... which has Pt>10. '
        #print '                  selectedPFTaus:  >= 1 PFTau which pass fixedConeHighEffPFTauDiscriminationByLeadingPionPtCut'
        print '                selectedCaloTaus:  >= 1 CaloTau which pass '
        print '                                   caloRecoTauDiscriminationByLeadingTrackPtCut (lead trk pt>5.)'
        print '                   selectedMuons:  >= 1 global muon with pt>10 & abs(eta)< 2.5 '
        #print '                    muPFTauPairs:  produce mu-PFTau pairs '
        #print '            selectedMuPFTauPairs:  require mu-PFTau pair separated by R>0.7 '
        print '                  muCaloTauPairs:  produce mu-CaloTau pairs '
        print '          selectedMuCaloTauPairs:  require mu-CaloTau pair separated by R>0.7 '
        print '-----------------------------------------------------------------------'
            
if __name__ == '__main__' :
    main(sys.argv[1:])

