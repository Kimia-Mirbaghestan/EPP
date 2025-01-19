import ROOT
import ctypes
import copy
import math
from Samples import samp

class MyAnalysis(object):
   
    def __init__(self, sample):

        """ The Init() function is called when an object MyAnalysis is initialised
        The tree corresponding to the specific sample is picked up
        and histograms are booked.
        """
        self._tree = ROOT.TTree()
        self.EventPassed = 0
        self.Triggered = 0
        if(sample not in samp.keys() and sample != "data"):
            print ("Error") #RuntimeError("Sample %s not valid. please, choose among these: %s" % (sample, str(samp.keys())) )
            exit
        self.histograms = {}
        self.sample = sample
        self._file = ROOT.TFile("files/"+sample+".root")
        self._file.cd()
        tree = self._file.Get("events")
        self._tree = tree
        self.acceptance = 0
        self.integral_aftercut = 0
        self.integral_beforecut = 0
        self.triggered_pass = 0
        self.acceptance_uncertainty = 0
        self.Trigger_Efficiency = 1.0
        self.Trigger_Efficiency_uncertainty = 0
        self.nEvents = self._tree.GetEntries()
        print ("-------------------------------------------------------")
        print ("Number of entries for " + self.sample + ": " + str(self.nEvents))
        ### Book histograms
        self.bookHistos()

    def getTree(self):
        return self._tree

    def getHistos(self):
        return self.histograms

    def bookHistos(self):
        h_nJet = ROOT.TH1F("NJet","#of jets", 6, -0.5, 6.5)
        h_nJet.SetXTitle("%# of jets")
        self.histograms["NJet"] = h_nJet

        h_nJetFinal = ROOT.TH1F("NJetFinal","#of jets", 6, -0.5, 6.5) #Plotted
        h_nJetFinal.SetXTitle("%# of jets")
        self.histograms["NJetFinal"] = h_nJetFinal

        h_MuonIso = ROOT.TH1F("Muon_Iso","Muon Isolation", 25, 0., 3.)
        h_MuonIso.SetXTitle("Muon Isolation")
        self.histograms["Muon_Iso"] = h_MuonIso

        h_NIsoMu = ROOT.TH1F("NIsoMu","Number of isolated muons", 5, 0.5, 5.5) #Plotted
        h_NIsoMu.SetXTitle("Number of isolated muons")
        self.histograms["NIsoMu"] = h_NIsoMu

        h_MuonPt = ROOT.TH1F("Muon_Pt","Muon P_T", 50, 0., 200.) #Plotted
        h_MuonPt.SetXTitle("Muon P_T")
        self.histograms["Muon_Pt"] = h_MuonPt

        h_METpt = ROOT.TH1F("MET_Pt","MET P_T", 25, 0., 300.) #Plotted
        h_METpt.SetXTitle("MET P_T")
        self.histograms["MET_Pt"] = h_METpt

        h_JetPt = ROOT.TH1F("Jet_Pt","Jet P_T", 50, 0., 200.) #Plotted
        h_JetPt.SetXTitle("Jet P_T")
        self.histograms["Jet_Pt"] = h_JetPt

        h_JetBtag = ROOT.TH1F("Jet_Btag","Jet B tag", 10, 1., 6.)
        h_JetBtag.SetXTitle("Jet B tag")
        self.histograms["Jet_btag"] = h_JetBtag

        h_NBtag = ROOT.TH1F("NBtag","Jet B tag", 12, 0.5, 12.5) #Plotted
        h_NBtag.SetXTitle("Number of B tagged jets")
        self.histograms["NBtag"] = h_NBtag
        
        h_nJet_nocut = ROOT.TH1F("NJet_NoCut","# of Jet", 12, -0.5, 12.5) #To Check Integral without Cut
        h_nJet_nocut.SetXTitle("Number of Jets without any Cut")
        self.histograms["NJet_NoCut"] =  h_nJet_nocut
        
        h_nJet_triggered = ROOT.TH1F("Triggered","# of Jet", 12, -0.5, 12.5) #To Check Integral without Cut
        h_nJet_triggered.SetXTitle("Number of Jets with Trigger")
        self.histograms["Triggered"] =  h_nJet_triggered

    def saveHistos(self):
        outfilename = "Histogram_Root/"+self.sample + "_histos.root"
        outfile = ROOT.TFile(outfilename, "RECREATE")
        outfile.cd()
        for h in self.histograms.values():
            h.Write()
        outfile.Close()

    ### processEvent function implements the actions to perform on each event
    ### This is the place where to implement the analysis strategy: study of most sensitive variables
    ### and signal-like event selection
    
    def processEvent2(self, entry):
        tree = self.getTree()
        tree.GetEntry(entry)
        w = tree.EventWeight
        nJets = tree.NJet
        self.histograms["NJet_NoCut"].Fill(nJets,w)
        
        

    def processEvent(self, entry):
        tree = self.getTree()
        tree.GetEntry(entry)
        w = tree.EventWeight

        ### Muon selection - Select events with at least 1 isolated muon
        ### with pt>25 GeV to match trigger requirements
        muonPtCut = 25.
        muonRelIsoCut = 0.05
        nIsoMu = 0
        muon_eta_cut = 2.1
        
        
        
        # Object selection: selecting the muon
        for m in range(tree.NMuon):
            muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
            self.histograms["Muon_Iso"].Fill(tree.Muon_Iso[m], w)
            if muon.Pt() > muonPtCut and (tree.Muon_Iso[m]/muon.Pt()) < muonRelIsoCut:
                if (abs(muon.Eta()) < muon_eta_cut):
                   nIsoMu += 1
        
        # Jet and b-jet counts
        nJets = 0
        nBjets = 0
        metPx=tree.MET_px
        metPy=tree.MET_py
        metPt=math.sqrt(metPx**2 + metPy**2)
        
        # Object selection: counting jets and b-jets
        for j in range(tree.NJet):
            jetPtCut = 35.0  # Define a jet pt cut
            jetBtagCut = 2.0  # Define a b-tagging cut
            single_topCut = 2.0 #Define single top cut
        
            jet = ROOT.TLorentzVector(tree.Jet_Px[j], tree.Jet_Py[j], tree.Jet_Pz[j], tree.Jet_E[j])
            if jet.Pt() > jetPtCut:
                nJets += 1
                # Check if the jet is b-tagged
                if tree.Jet_btag[j] > jetBtagCut:
                    nBjets += 1
              
        # Event Selection (Apply cuts)
        if nJets >= 4 and nIsoMu >= 1 and nBjets >= 2 and metPt > 30:
            self.EventPassed += 1
            if tree.triggerIsoMu24 :
                self.Triggered += 1
                self.histograms["Triggered"].Fill(nJets,w)
            self.histograms["NJetFinal"].Fill(nJets, w)
            self.histograms["Muon_Pt"].Fill(muon.Pt(), w)
            self.histograms["NIsoMu"].Fill(nIsoMu, w)
            self.histograms["MET_Pt"].Fill(metPt,w)
            self.histograms["Jet_Pt"].Fill(jet.Pt(),w)
            self.histograms["NBtag"].Fill(nBjets,w)
            

    def processEvents(self):
        minBin_beforecut = self.histograms["NJet_NoCut"].GetXaxis().GetFirst()
        maxBin_beforecut = self.histograms["NJet_NoCut"].GetXaxis().GetLast()
        nDataErr_beforecut = ctypes.c_double(0.0)
        minBin_aftercut = self.histograms["NBtag"].GetXaxis().GetFirst()
        maxBin_aftercut = self.histograms["NBtag"].GetXaxis().GetLast()
        nDataErr_aftercut = ctypes.c_double(0.0)
        minBin_triggered = self.histograms["Triggered"].GetXaxis().GetFirst()
        maxBin_triggered = self.histograms["Triggered"].GetXaxis().GetLast()
        nDataErr_triggered = ctypes.c_double(0.0)
        
        nevts = self.nEvents
        for i in range(nevts):
            self.processEvent(i)
            if True:
                self.processEvent2(i)
        if True:
            self.integral_beforecut = self.histograms["NJet_NoCut"].IntegralAndError(minBin_beforecut,maxBin_beforecut,nDataErr_beforecut)
            self.integral_aftercut = self.histograms["NBtag"].IntegralAndError(minBin_aftercut,maxBin_aftercut,nDataErr_aftercut)
            self.acceptance = (self.integral_aftercut/self.integral_beforecut) 
            self.acceptance_uncertainty = (nDataErr_aftercut.value/self.integral_beforecut)**2 +((nDataErr_beforecut.value*self.integral_aftercut)/(self.integral_beforecut**2))**2
            print ("Integral without Cut : "+str(self.integral_beforecut)+" Integral uncertainty without Cut : "+str(nDataErr_beforecut.value))
            print ("Integral with Cut : "+str(self.integral_aftercut)+" Integral uncertainty with Cut : "+str(nDataErr_aftercut.value))
            print ("Acceptance : "+str(self.acceptance)+" acceptance uncertanty : " +str(self.acceptance_uncertainty))


        print ("Events Passed : " + str(self.EventPassed))
        if self.sample == "ttbar":
            self.triggered_pass = self.histograms["Triggered"].IntegralAndError(minBin_triggered,maxBin_triggered,nDataErr_triggered)
            print ("Events Passed with Trigger : " + str(self.triggered_pass))
            self.Trigger_Efficiency = (self.triggered_pass/self.integral_aftercut) 
            self.Trigger_Efficiency_uncertainty = (nDataErr_triggered.value/self.integral_aftercut)**2 +((nDataErr_triggered.value*self.triggered_pass)/(self.integral_aftercut**2))**2
            print ("Trigger Efficiency : " + str(self.Trigger_Efficiency)+" Trigger Efficiency uncertainty : " + str(self.Trigger_Efficiency_uncertainty))
            self.cross_section()
        print ("Cut Efficiency : " + str((self.EventPassed/nevts)*100))
        self.saveHistos()
        
        
        
        
    def getNTriggers(self):
        return self.n_triggers
        

    
    
    def cross_section(self):
        # Assuming the number of observed events is stored in the histogram "NJetFinal"
        N_observed = 39    #45           #45 * (36941/214)  #7767.96

        # Integrated luminosity (in pb^-1)
        luminosity = 50  # Example value, replace with the actual integrated luminosity
        

        # Compute cross-section
        
        sigma_observed = N_observed / (luminosity * self.acceptance* self.Trigger_Efficiency)
        
        un_N_observed = math.sqrt(N_observed)
        un_luminosity = 0.1
        un_Trigger_Efficiency = self.Trigger_Efficiency_uncertainty
        un_acceptance = self.acceptance_uncertainty
        
        total_uncertainty = math.sqrt((un_N_observed/luminosity * self.acceptance* self.Trigger_Efficiency)**2 + ((un_luminosity*N_observed)/((luminosity**2)*self.acceptance*self.Trigger_Efficiency))**2 + ((un_Trigger_Efficiency*N_observed)/((self.acceptance**2)*self.Trigger_Efficiency*luminosity*luminosity))**2 + ((un_acceptance*N_observed)/((self.Trigger_Efficiency)*self.acceptance*luminosity))**2)
        
        sys_un = un_N_observed/luminosity * self.acceptance* self.Trigger_Efficiency
        stat_un = math.sqrt(((un_luminosity*N_observed)/((luminosity**2)*self.acceptance*self.Trigger_Efficiency))**2 + ((un_Trigger_Efficiency*N_observed)/((self.acceptance**2)*self.Trigger_Efficiency*luminosity*luminosity))**2 + ((un_acceptance*N_observed)/((self.Trigger_Efficiency)*self.acceptance*luminosity))**2)
        
        print("Observed cross-section:", sigma_observed, "pb")
        print("total uncertainty:", total_uncertainty, "pb")
        print("sys uncertainty:", sys_un, "pb")
        print("stat uncertainty:", stat_un, "pb")







