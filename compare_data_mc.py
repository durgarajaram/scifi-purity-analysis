import ROOT
import sys

def compare_plots(fmc, fdata):

    _mc = ROOT.TFile(fmc)
    _dt = ROOT.TFile(fdata)

    kl = [ key.GetName() for key in ROOT.gDirectory.GetListOfKeys() ]
    print 'list of keys: ',kl

    c = ROOT.TCanvas("c", "c", 800, 800)
    c.SetGridx()
    c.SetGridy()

    ROOT.gStyle.SetPadGridX(1)
    ROOT.gStyle.SetPadGridY(1)

    ####################################################

    mcnpe = _mc.Get("upstream/cluster_npe")
    mcnpe.SetLineColor(2)
    dcnpe = _dt.Get("upstream/cluster_npe")

    nd = dcnpe.GetEntries()
    nm = mcnpe.GetEntries()
    
    dcnpe.Scale(nm/nd)

    dmax = dcnpe.GetMaximum()
    mmax = mcnpe.GetMaximum()

    max = dcnpe.GetMaximum()*1.05 if dcnpe.GetMaximum() > mcnpe.GetMaximum() else mcnpe.GetMaximum()*1.05

    dcnpe.SetMaximum(max)
    mcnpe.SetMaximum(max)

    leg = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    leg.SetFillColor(0)
    leg.AddEntry(dcnpe, "Data")
    leg.AddEntry(mcnpe, "MC")
    dcnpe.Draw()
    mcnpe.Draw("esame")
    leg.Draw()
    c.Print("upstream_cluster_npe.png")

    ####################################################

    mcnpe = _mc.Get("upstream/spacepoint_npe")
    mcnpe.SetLineColor(2)
    dcnpe = _dt.Get("upstream/spacepoint_npe")

    nd = dcnpe.GetEntries()
    nm = mcnpe.GetEntries()
    
    dcnpe.Scale(nm/nd)

    dmax = dcnpe.GetMaximum()
    mmax = mcnpe.GetMaximum()

    max = dcnpe.GetMaximum()*1.05 if dcnpe.GetMaximum() > mcnpe.GetMaximum() else mcnpe.GetMaximum()*1.05

    dcnpe.SetMaximum(max)
    mcnpe.SetMaximum(max)

    leg = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    leg.SetFillColor(0)
    leg.AddEntry(dcnpe, "Data")
    leg.AddEntry(mcnpe, "MC")
    dcnpe.Draw()
    mcnpe.Draw("esame")
    leg.Draw()
    c.Print("upstream_spacepoint_npe.png")

    ####################################################
    mcnpe = _mc.Get("downstream/cluster_npe")
    mcnpe.SetLineColor(2)
    dcnpe = _dt.Get("downstream/cluster_npe")

    nd = dcnpe.GetEntries()
    nm = mcnpe.GetEntries()
    
    dcnpe.Scale(nm/nd)

    max = dcnpe.GetMaximum()*1.05 if dcnpe.GetMaximum() > mcnpe.GetMaximum() else mcnpe.GetMaximum()*1.05

    dcnpe.SetMaximum(max)
    mcnpe.SetMaximum(max)

    dcnpe.Draw()
    mcnpe.Draw("esame")
    leg.Draw()
    c.Print("downstream_cluster_npe.png")

    ####################################################
    mcnpe = _mc.Get("downstream/spacepoint_npe")
    mcnpe.SetLineColor(2)
    dcnpe = _dt.Get("downstream/spacepoint_npe")

    nd = dcnpe.GetEntries()
    nm = mcnpe.GetEntries()
    
    dcnpe.Scale(nm/nd)

    max = dcnpe.GetMaximum()*1.05 if dcnpe.GetMaximum() > mcnpe.GetMaximum() else mcnpe.GetMaximum()*1.05

    dcnpe.SetMaximum(max)
    mcnpe.SetMaximum(max)

    dcnpe.Draw()
    mcnpe.Draw("esame")
    leg.Draw()
    c.Print("downstream_spacepoint_npe.png")

    for s in range(5):
        ################################################
        c.Clear()
        name = "sp_npe_used_"+str(s+1)
        of = "upstream_"+name+".png"
#       mhname = "h"+name
        name = "upstream/"+name
        print 'Getting..',name

        mhname = _mc.Get(name)
        mhname.SetLineColor(2)
        
        dhname = _dt.Get(name)

        nd = dhname.GetEntries()
        nm = mhname.GetEntries()
        dhname.Scale(nm/nd)

        max = dhname.GetMaximum()*1.05 if dcnpe.GetMaximum() > mhname.GetMaximum() else mcnpe.GetMaximum()*1.05

        dhname.SetMaximum(max)
        mhname.SetMaximum(max)

        mhname.Draw()
        dhname.Draw("esame")
        c.Print(of)

        ################################################
        c.Clear()
        name = "sp_npe_used_"+str(s+1)
        of = "downstream_"+name+".png"
#       mhname = "h"+name
        name = "downstream/"+name
        print 'Getting..',name

        mhname = _mc.Get(name)
        mhname.SetLineColor(2)
        
        dhname = _dt.Get(name)

        nd = dhname.GetEntries()
        nm = mhname.GetEntries()
        dhname.Scale(nm/nd)


        max = dhname.GetMaximum()*1.05 if dcnpe.GetMaximum() > mhname.GetMaximum() else mcnpe.GetMaximum()*1.05

        dhname.SetMaximum(max)
        mhname.SetMaximum(max)

        mhname.Draw()
        dhname.Draw("esame")
        c.Print(of)

        ###############################################
        c.Clear()
        c.Divide(1,2)
            
        name = "sp_xy_used_"+str(s+1)
        of = "upstream_"+name+".png"
        name = "upstream/"+name

        mhname = _mc.Get(name)
        mhname.SetLineColor(2)
        
        dhname = _dt.Get(name)

        nd = dhname.GetEntries()
        nm = mhname.GetEntries()

        dhname.Scale(nm/nd)

        c.cd(1)
        mhname.SetTitle("SP used MC")
        mhname.Draw("colz")
        c.cd(2)
        dhname.SetTitle("SP used Data")
        dhname.Draw("colz")

        c.Print(of)
        ###############################################
        c.Clear()
        c.Divide(1,2)
            
        name = "sp_xy_used_"+str(s+1)
        of = "downstream_"+name+".png"
        name = "downstream/"+name

        mhname = _mc.Get(name)
        mhname.SetLineColor(2)
        
        dhname = _dt.Get(name)

        nd = dhname.GetEntries()
        nm = mhname.GetEntries()

        dhname.Scale(nm/nd)

        c.cd(1)
        mhname.SetTitle("SP used MC")
        mhname.Draw("colz")
        c.cd(2)
        dhname.SetTitle("SP used Data")
        dhname.Draw("colz")

        c.Print(of)

        for p in range(3):
            ################################################
            c.Clear()
            name = "digit_npe_s"+str(s+1)+"_p"+str(p)
            of = "upstream_"+name+".png"
            name = "upstream/"+name
            print 'Getting..',name

            mhname = _mc.Get(name)
            mhname.SetLineColor(2)
        
            dhname = _dt.Get(name)

            nd = dhname.GetEntries()
            nm = mhname.GetEntries()

            dhname.Scale(nm/nd)
  
            max = dhname.GetMaximum()*1.05 if dcnpe.GetMaximum() > mhname.GetMaximum() else mcnpe.GetMaximum()*1.05
            print dhname.GetMaximum(), mhname.GetMaximum()
            print max

            dhname.SetMaximum(max)
            mhname.SetMaximum(max)

            mhname.Draw()
            dhname.Draw("esame")
            c.Print(of)
            ##############################################
            c.Clear()
            name = "digit_npe_s"+str(s+1)+"_p"+str(p)
            of = "downstream_"+name+".png"
            name = "downstream/"+name
            print 'Getting..',name

            mhname = _mc.Get(name)
            mhname.SetLineColor(2)
        
            dhname = _dt.Get(name)

            nd = dhname.GetEntries()
            nm = mhname.GetEntries()

            dhname.Scale(nm/nd)
            max = dhname.GetMaximum()*1.05 if dcnpe.GetMaximum() > mhname.GetMaximum() else mcnpe.GetMaximum()*1.05

            print dhname.GetMaximum(), mhname.GetMaximum()
            print max
            dhname.SetMaximum(max)
            mhname.SetMaximum(max)

            mhname.Draw()
            dhname.Draw("esame")
            c.Print(of)
            ##############################################
            c.Clear()
            c.Divide(1,2)

            name = "digit_npeVchan_s"+str(s+1)+"_p"+str(p)
            of = "upstream"+name+".png"
            name = "upstream/"+name
            print 'Getting..',name

            mhname = _mc.Get(name)
            mhname.SetLineColor(2)
        
            dhname = _dt.Get(name)

            nd = dhname.GetEntries()
            nm = mhname.GetEntries()

            c.cd(1)
            mhname.SetTitle("npe vs chan MC")
            mhname.Draw("colz")
            c.cd(2)
            dhname.SetTitle("npe vs chan Data")
            dhname.Draw("colz")
            c.Print(of)
            ##############################################
            c.Clear()
            c.Divide(1,2)

            name = "digit_npeVchan_s"+str(s+1)+"_p"+str(p)
            of = "downstream"+name+".png"
            name = "downstream/"+name
            print 'Getting..',name

            mhname = _mc.Get(name)
            mhname.SetLineColor(2)
        
            dhname = _dt.Get(name)

            nd = dhname.GetEntries()
            nm = mhname.GetEntries()

            c.cd(1)
            mhname.SetTitle("npe vs chan MC")
            mhname.Draw("colz")
            c.cd(2)
            dhname.SetTitle("npe vs chan Data")
            dhname.Draw("colz")
            c.Print(of)
            ##############################################
            c.Clear()
            name = "cluster_npe_s"+str(s+1)+"_p"+str(p)
            of = "upstream_"+name+".png"
            name = "upstream/"+name
            print 'Getting..',name

            mhname = _mc.Get(name)
            mhname.SetLineColor(2)
        
            dhname = _dt.Get(name)

            nd = dhname.GetEntries()
            nm = mhname.GetEntries()

            dhname.Scale(nm/nd)
  
            max = dhname.GetMaximum()*1.05 if dcnpe.GetMaximum() > mhname.GetMaximum() else mcnpe.GetMaximum()*1.05

            dhname.SetMaximum(max)
            mhname.SetMaximum(max)

            mhname.Draw()
            dhname.Draw("esame")
            c.Print(of)
            ##############################################
            c.Clear()
            name = "cluster_npe_s"+str(s+1)+"_p"+str(p)
            of = "downstream_"+name+".png"
            name = "downstream/"+name
            print 'Getting..',name

            mhname = _mc.Get(name)
            mhname.SetLineColor(2)
        
            dhname = _dt.Get(name)

            nd = dhname.GetEntries()
            nm = mhname.GetEntries()

            dhname.Scale(nm/nd)
  
            max = dhname.GetMaximum()*1.05 if dcnpe.GetMaximum() > mhname.GetMaximum() else mcnpe.GetMaximum()*1.05

            dhname.SetMaximum(max)
            mhname.SetMaximum(max)

            mhname.Draw()
            dhname.Draw("esame")
            c.Print(of)
            ##############################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print 'Usage: ',__file__,'mc-file-name data-file-name'
        exit(1)
    mc_filename = sys.argv[1]
    data_filename = sys.argv[2]
    ROOT.gROOT.SetBatch(1)
    c = ROOT.TCanvas("c", "c", 800, 800)
    compare_plots(mc_filename, data_filename)

