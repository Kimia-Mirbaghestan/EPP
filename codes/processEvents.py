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
